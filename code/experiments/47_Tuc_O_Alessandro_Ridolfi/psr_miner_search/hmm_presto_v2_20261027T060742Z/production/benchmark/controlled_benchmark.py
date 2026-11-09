#!/usr/bin/env python3
"""Measure saved lean HMM configurations and a complete PRESTO baseline."""
import argparse
import ast
import csv
from decimal import Decimal
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import statistics
import sys
import time
import re

from presto_pilot import DATA_ROOT, EXPERIMENT, inventory, product, read_info, timed


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def load_params(path):
    values = {}
    for line in path.read_text().splitlines():
        if "=" in line and not line.lstrip().startswith("#"):
            key, value = line.split("=", 1)
            values[key.strip()] = ast.literal_eval(value.strip())
    return values


def instrument(source):
    # Insert observations only. Preserve every original numerical calculation.
    preamble = '''import time as _bench_time, json as _bench_json, os as _bench_os
def _bench_mark(name, **details):
    row = dict(event=name, wall=_bench_time.perf_counter(), cpu=_bench_time.process_time(), **details)
    with open(_bench_os.environ["HMM_BENCHMARK_PHASE_FILE"], "a") as handle:
        handle.write(_bench_json.dumps(row) + "\\n")
_bench_mark("script_start")
'''
    insertions = [
        ("    data = np.fromfile(args.infile, dtype=np.single)", '    _bench_mark("read_start")\n', "before"),
        ("    Tfull = len(data) * tsamp", '    _bench_mark("read_done", real_samples=len(data))\n', "before"),
        ("    _, fs = compute_sft(downsample_data, downsample_tsamp, Tsft,", '    _bench_mark("spectrogram_start")\n', "before"),
        ("    delta, backptrs = viterbi(spec)", '    _bench_mark("spectrogram_done", segments=Nsft, Tsft=Tsft, frequency_bins=len(fs))\n', "before"),
        ("    delta, backptrs = viterbi(spec)", '\n    _bench_mark("tracking_done")', "after"),
        ("    if args.lean:", '    _bench_mark("products_done")\n', "before"),
    ]
    for target, addition, side in insertions:
        if source.count(target) != 1:
            raise ValueError("Source does not match the reviewed instrumentation point: " + target)
        source = source.replace(target, addition + target if side == "before" else target + addition, 1)
    result = preamble + source
    ast.parse(result)
    return result


def machine():
    return {"hostname": os.uname().nodename, "python": sys.version,
            "affinity": sorted(os.sched_getaffinity(0)),
            "cpuinfo": next((s for s in Path("/proc/cpuinfo").read_text().splitlines() if s.startswith("model name")), ""),
            "job": os.environ.get("SLURM_JOB_ID"),
            "versions": {name: importlib.metadata.version(name) for name in ["numpy","scipy","matplotlib"] if name in {d.metadata["Name"].lower() for d in importlib.metadata.distributions()} }}


def run_hmm(args):
    directory = args.root / "hmm"
    directory.mkdir(parents=True, exist_ok=False)
    sources = Path(__file__).resolve().parent / "hmm_sources"
    original = sources / "viterbi_pipeline.py"
    instrumented = directory / "instrumented_pipeline.py"
    instrumented.write_text(instrument(original.read_text()))
    shutil.copyfile(sources / "viterbi.py", directory / "viterbi.py")
    save_json(directory / "machine.json", machine())
    save_json(directory / "source_hashes.json", {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sources.glob("*.py")})
    fields = ["label_Nt", "repeat", "actual_segments", "Tsft_s", "covered_s", "command_wall_s", "command_cpu_s", "read_wall_s", "spectrogram_wall_s", "tracking_wall_s", "products_wall_s", "measured_peak_L"]
    measurements = []
    with (directory / "hmm_timings.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for repeat in range(1, args.repeats + 1):
            for nt in [16, 32, 64, 128, 256, 512, 1024]:
                run = directory / ("Nt{}_repeat{}".format(nt, repeat))
                run.mkdir()
                source=DATA_ROOT/"DM24.40_bary/47Tuc_DM24.40_bary.dat"
                info=read_info(source.with_suffix(".inf"))
                values=dict(infile=str(source), tsamp=float(info["dt_s"]), segments=nt,
                    f0=374.0,bw=10.0,padding_factor=1,num_harm=1,top_paths=1,save_delta=False,
                    dm=24.4,mjd_start=float(info["epoch_mjd"]))
                values["out_prefix"] = str(run / "search")
                params = run / "search.params"
                params.write_text("\n".join("{} = {!r}".format(k, v) for k, v in values.items()) + "\n")
                os.environ["HMM_BENCHMARK_PHASE_FILE"] = str(run / "phases.jsonl")
                result = timed("hmm", [sys.executable, str(instrumented), "--params", str(params), "--lean"], run)
                phases = {x["event"]: x for x in [json.loads(line) for line in (run / "phases.jsonl").read_text().splitlines()]}
                def elapsed(start, end):
                    return phases[end]["wall"] - phases[start]["wall"]
                actual = phases["spectrogram_done"]
                def peak(path):
                    return max(float(line.split()[1]) for line in path.read_text().splitlines() if line.strip() and not line.startswith("#"))
                sampling=json.loads((run/"search_sampling.json").read_text())
                if sampling["actual_segments"]!=nt:
                    raise ValueError("Incorrect segment count")
                if sampling["frequency_axis_version"]!=2 or sampling["unused_trailing_samples"]>=nt:
                    raise ValueError("Invalid corrected sampling metadata")
                old_peak = 0.0
                new_peak = peak(run / "search_loglike_curve.dat")
                row = {"label_Nt": nt, "repeat": repeat, "actual_segments": actual["segments"], "Tsft_s": actual["Tsft"],
                       "covered_s": actual["segments"] * actual["Tsft"], "command_wall_s": result["wall_s"], "command_cpu_s": result["cpu_s"],
                       "read_wall_s": elapsed("read_start", "read_done"), "spectrogram_wall_s": elapsed("spectrogram_start", "spectrogram_done"),
                       "tracking_wall_s": elapsed("spectrogram_done", "tracking_done"), "products_wall_s": elapsed("tracking_done", "products_done"),
                       "measured_peak_L": new_peak}
                writer.writerow(row)
                handle.flush()
                save_json(run / "measurement.json", row)
                measurements.append(row)
                print("HMM_RESULT " + json.dumps(row), flush=True)
    save_json(directory / "summary.json", {"status": "complete", "runs": len(measurements),
        "corrected_segments_verified": all(r["actual_segments"] == r["label_Nt"] for r in measurements),
        "by_label_Nt": {str(nt): {"median_command_wall_s": statistics.median(r["command_wall_s"] for r in measurements if r["label_Nt"] == nt),
                                     "median_command_cpu_s": statistics.median(r["command_cpu_s"] for r in measurements if r["label_Nt"] == nt)} for nt in [16,32,64,128,256,512,1024]}})


def extract_section(args):
    import numpy as np
    source = args.source
    info = read_info(source.with_suffix(".inf"))
    if args.parts < 1 or not 0 <= args.index < args.parts or source.stat().st_size != info["nsamples"] * 4:
        raise ValueError("Invalid partition index or source byte count")
    valid = info["onoff"][0][1] + 1 if info["breaks"] else info["nsamples"]
    if info["breaks"] and info["onoff"][0][0] != 0:
        raise ValueError("Input does not start with the expected continuous valid interval")
    if info["breaks"] and info["onoff"] != [[0, valid-1],[info["nsamples"]-1,info["nsamples"]-1]]:
        raise ValueError("Unexpected input gaps; refuse to reinterpret as trailing padding")
    first = valid * args.index // args.parts
    last = valid * (args.index + 1) // args.parts
    count = last - first
    if info["nsamples"] % args.parts:
        raise ValueError("Stored FFT length is not divisible by requested partition count")
    fft_count = info["nsamples"] // args.parts
    if fft_count < count or info["barycentric"] != 1:
        raise ValueError("Invalid FFT length or barycentric flag")
    before = source.stat()
    inp = np.memmap(str(source), dtype=np.float32, mode="r", shape=(info["nsamples"],))
    mean = float(np.mean(inp[first:last], dtype=np.float64))
    dst = product(args.destination, ".dat")
    with dst.open("xb"):
        pass
    out = np.memmap(str(dst), dtype=np.float32, mode="r+", shape=(fft_count,)) if dst.stat().st_size else None
    if out is None:
        # Allocate only a new, exclusively created working file.
        with dst.open("r+b") as f:
            f.truncate(fft_count * 4)
        out = np.memmap(str(dst), dtype=np.float32, mode="r+", shape=(fft_count,))
    block = 1024 * 1024
    for offset in range(0, count, block):
        length = min(block, count - offset)
        np.subtract(inp[first + offset:first + offset + length], np.float32(mean), out=out[offset:offset + length])
    out[count:] = 0
    out.flush()
    del out, inp
    lines = []
    for line in source.with_suffix(".inf").read_text().splitlines():
        s = line.lstrip()
        if s.startswith("On/Off bin pair"):
            continue
        if s.startswith("Data file name without suffix"):
            line = " Data file name without suffix          =  " + args.destination.name
        elif s.startswith("Epoch of observation"):
            epoch = Decimal(info["epoch_mjd"]) + Decimal(first) * Decimal(info["dt_s"]) / Decimal(86400)
            line = " Epoch of observation (MJD)             =  " + str(epoch)
        elif s.startswith("Number of bins in the time series"):
            line = " Number of bins in the time series      =  " + str(fft_count)
        elif s.startswith("Any breaks in the data?"):
            line = " Any breaks in the data? (1 yes, 0 no)  =  " + str(int(fft_count > count))
            if fft_count > count:
                line += "\n On/Off bin pair #  1                   =  0 , {}\n On/Off bin pair #  2                   =  {} , {}".format(count - 1, fft_count - 1, fft_count - 1)
        lines.append(line)
    product(args.destination, ".inf").write_text("\n".join(lines) + "\n")
    after = source.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise RuntimeError("Source changed while preparing section")
    save_json(product(args.destination, ".section.json"), {
        "source": str(source), "dm": info["dm"], "parts": args.parts, "index": args.index,
        "start_sample": first, "end_sample_exclusive": last, "valid_samples": count,
        "fft_samples": fft_count, "valid_duration_s": float(Decimal(count) * Decimal(info["dt_s"])),
        "fft_duration_s": float(Decimal(fft_count) * Decimal(info["dt_s"])),
        "subtracted_mean_float32": float(np.float32(mean)), "padding": "zeros after mean subtraction"})


def parse_presto_candidates(path, dm, parts, index):
    rows=[]
    for line in path.read_text().splitlines():
        if "Power /" in line: break  # The next table describes individual harmonics.
        cols=line.split()
        if len(cols)<11 or not cols[0].isdigit(): continue
        def numeric(text):
            cleaned=re.sub(r"\([^)]*\)", "",text).replace("x10^","e")
            return float(cleaned)
        rows.append(dict(candidate_id=f"presto_DM{dm:.2f}_parts{parts}_section{index:02d}_cand{int(cols[0]):06d}",
            dm_best=dm,parts=parts,section=index,presto_candidate_number=int(cols[0]),
            sigma=numeric(cols[1]),summed_power=numeric(cols[2]),coherent_power=numeric(cols[3]),
            numharm=int(cols[4]),peak_freq_hz=numeric(cols[6]),fdot_hz_s=numeric(cols[8]),
            z=numeric(cols[9]),source_table=str(path),notes=" ".join(cols[11:]),retained=True,rfi_policy="flag_only"))
    return rows


def run_presto(args):
    directory = args.root / "presto"
    directory.mkdir(parents=True, exist_ok=False)
    programs = {name: shutil.which(name) for name in ["realfft", "rednoise", "accelsearch"]}
    if not all(programs.values()):
        raise ValueError("A required PRESTO executable is missing")
    trials = inventory(DATA_ROOT)
    save_json(directory / "input_manifest.json", list(trials.values()))
    save_json(directory / "machine.json", dict(machine(), programs=programs, fhi=1059, numharm=1, zmax=200, wmax=0))
    selected = sorted(trials) if args.all_dms else [2440]
    measurements = []
    candidate_rows = []
    fields = ["dm", "parts", "section", "valid_duration_s", "fft_duration_s", "preparation_wall_s", "preparation_cpu_s", "search_wall_s", "search_cpu_s", "fft_wall_s", "whitening_wall_s", "acceleration_wall_s", "total_wall_s", "total_cpu_s", "max_rss_kib"]
    with (directory / "presto_timings.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for dm_key in selected:
            source = Path(trials[dm_key]["dat"])
            for parts in [1, 8, 16]:
                for index in range(parts):
                    run = directory / "DM{:.2f}".format(dm_key / 100) / "parts{}".format(parts) / "section{:02d}".format(index)
                    run.mkdir(parents=True)
                    prefix = run / "section"
                    commands = [
                        ("01_preparation", [sys.executable, str(Path(__file__).resolve()), "extract", str(source), str(prefix), "--parts", str(parts), "--index", str(index)]),
                        ("02_fft", [programs["realfft"], "-fwd", "section.dat"]),
                        ("03_whitening", [programs["rednoise"], "section.fft"]),
                        ("04_acceleration", [programs["accelsearch"], "-ncpus", "1", "-numharm", "1", "-zmax", "200", "-wmax", "0", "-sigma", "2", "-flo", "50", "-fhi", "1059", "section_red.fft"]),
                    ]
                    stages = [timed(label, command, run) for label, command in commands]
                    candidate_table=run/"section_red_ACCEL_200"
                    if not candidate_table.is_file(): raise FileNotFoundError(candidate_table)
                    candidates=parse_presto_candidates(candidate_table,dm_key/100,parts,index)
                    candidate_rows.extend(candidates)
                    save_json(run/"candidates.json",candidates)
                    section = json.loads(product(prefix, ".section.json").read_text())
                    row = {"dm": dm_key / 100, "parts": parts, "section": index,
                           "valid_duration_s": section["valid_duration_s"], "fft_duration_s": section["fft_duration_s"],
                           "preparation_wall_s": stages[0]["wall_s"], "preparation_cpu_s": stages[0]["cpu_s"],
                           "search_wall_s": sum(x["wall_s"] for x in stages[1:]), "search_cpu_s": sum(x["cpu_s"] for x in stages[1:]),
                           "fft_wall_s": stages[1]["wall_s"], "whitening_wall_s": stages[2]["wall_s"], "acceleration_wall_s": stages[3]["wall_s"],
                           "total_wall_s":sum(x["wall_s"] for x in stages),
                           "total_cpu_s":sum(x["cpu_s"] for x in stages),
                           "max_rss_kib": max(x["max_rss_kib"] for x in stages)}
                    save_json(run / "measurement.json", row)
                    measurements.append(row)
                    writer.writerow(row)
                    handle.flush()
                    print("PRESTO_RESULT " + json.dumps(row), flush=True)
                    # Remove only transient files created by this completed section.
                    # Retain metadata, candidates, logs, timings and the reconstruction recipe.
                    for name in ["section.dat", "section.fft", "section_red.fft"]:
                        (run / name).unlink()
    columns=["candidate_id","dm_best","parts","section","presto_candidate_number","sigma","summed_power","coherent_power","numharm","peak_freq_hz","fdot_hz_s","z","source_table","notes","retained","rfi_policy"]
    with (directory/"candidates_all.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=columns);writer.writeheader();writer.writerows(candidate_rows)
    if len({r["candidate_id"] for r in candidate_rows})!=len(candidate_rows): raise ValueError("Duplicate PRESTO candidate IDs")
    save_json(directory / "summary.json", {"status": "complete", "retained_candidates":len(candidate_rows), "candidate_threshold_note":"PRESTO sigma=2, not a calibrated match to HMM false alarm probability", "dm_trials": len(selected), "sections": len(measurements),
        "frequency_hz": [50, 1059], "numharm": 1, "zmax": 200, "wmax": 0,
        "by_parts": {str(parts): {"sections": sum(r["parts"] == parts for r in measurements),
                                 "preparation_wall_s": sum(r["preparation_wall_s"] for r in measurements if r["parts"] == parts),
                                 "preparation_cpu_s": sum(r["preparation_cpu_s"] for r in measurements if r["parts"] == parts),
                                 "search_wall_s": sum(r["search_wall_s"] for r in measurements if r["parts"] == parts),
                                 "search_cpu_s": sum(r["search_cpu_s"] for r in measurements if r["parts"] == parts)} for parts in [1,8,16]}})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="mode", required=True)
    hmm = sub.add_parser("hmm")
    hmm.add_argument("--root", type=Path, required=True)
    hmm.add_argument("--repeats", type=int, default=2)
    presto = sub.add_parser("presto")
    presto.add_argument("--root", type=Path, required=True)
    presto.add_argument("--all-dms", action="store_true")
    extract = sub.add_parser("extract")
    extract.add_argument("source", type=Path)
    extract.add_argument("destination", type=Path)
    extract.add_argument("--parts", type=int, required=True)
    extract.add_argument("--index", type=int, required=True)
    args = parser.parse_args()
    if args.mode != "extract" and not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Submit the controlled benchmark through Slurm")
    if args.mode == "hmm":
        run_hmm(args)
    elif args.mode == "presto":
        run_presto(args)
    else:
        extract_section(args)
