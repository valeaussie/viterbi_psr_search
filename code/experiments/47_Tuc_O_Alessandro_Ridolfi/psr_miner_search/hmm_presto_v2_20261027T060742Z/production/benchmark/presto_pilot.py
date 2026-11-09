#!/usr/bin/env python3
"""Time a single-core PRESTO pilot on an existing barycentred DM trial.

Uses Python's wait4 resource accounting, so GNU time is not required.
The extract subcommand only copies samples and adjusts their metadata.
"""

import argparse
import csv
import datetime
from decimal import Decimal
import hashlib
import json
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import sys
import time


EXPERIMENT = Path("/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi")
DATA_ROOT = EXPERIMENT / "47Tuc_blind_search_v1/stage2_dedisp"
WORK_ROOT = EXPERIMENT / "psr_miner_search"
FIELDS = {
    "Data file name without suffix": "basename",
    "Epoch of observation": "epoch_mjd",
    "Barycentered?": "barycentric",
    "Number of bins in the time series": "nsamples",
    "Width of each time series bin": "dt_s",
    "Any breaks in the data?": "breaks",
    "Dispersion measure": "dm",
}


def product(root, suffix):
    # A DM such as 24.40 is part of the basename, not a file extension.
    return Path(str(root) + suffix)


def read_info(path):
    info = {"onoff": []}
    for line in path.read_text().splitlines():
        text = line.lstrip()
        for prefix, name in FIELDS.items():
            if text.startswith(prefix):
                info[name] = line.rsplit("=", 1)[1].strip()
                break
        if text.startswith("On/Off bin pair"):
            info["onoff"].append([int(x.strip()) for x in line.rsplit("=", 1)[1].split(",")])
    for key in ("nsamples", "barycentric", "breaks"):
        info[key] = int(info[key])
    info["dm"] = float(info["dm"])
    return info


def fft_length_at_most(limit):
    """Largest 2/3/5/7-smooth length <= limit, divisible by 16."""
    values = {16}
    for factor in (2, 3, 5, 7):
        expanded = set()
        for value in values:
            while value <= limit:
                expanded.add(value)
                value *= factor
        values = expanded
    if not values:
        raise ValueError("Requested section is too short")
    return max(values)


def extract(source, destination, seconds):
    info = read_info(source.with_suffix(".inf"))
    if info["barycentric"] != 1:
        raise ValueError("Input is not marked barycentred")
    if source.stat().st_size != 4 * info["nsamples"]:
        raise ValueError("Input size does not match float32 samples in its .inf")
    count = fft_length_at_most(int(Decimal(str(seconds)) / Decimal(info["dt_s"])))
    if count > info["nsamples"]:
        raise ValueError("Requested section exceeds input duration")
    if info["breaks"] and not any(start == 0 and end >= count - 1 for start, end in info["onoff"]):
        raise ValueError("Requested first section includes a gap or trailing padding")
    digest = hashlib.sha256()
    remaining = 4 * count
    with source.open("rb") as src, product(destination, ".dat").open("xb") as dst:
        while remaining:
            block = src.read(min(4 * 1024 * 1024, remaining))
            if not block:
                raise ValueError("Unexpected end of input time series")
            dst.write(block)
            digest.update(block)
            remaining -= len(block)
    lines = []
    for line in source.with_suffix(".inf").read_text().splitlines():
        text = line.lstrip()
        if text.startswith("On/Off bin pair"):
            continue
        if text.startswith("Data file name without suffix"):
            line = " Data file name without suffix          =  " + destination.name
        elif text.startswith("Number of bins in the time series"):
            line = " Number of bins in the time series      =  " + str(count)
        elif text.startswith("Any breaks in the data?"):
            line = " Any breaks in the data? (1 yes, 0 no)  =  0"
        lines.append(line)
    product(destination, ".inf").write_text("\n".join(lines) + "\n")
    details = {
        "source": str(source), "start_sample": 0, "nsamples": count,
        "dt_s": info["dt_s"], "duration_s": float(count * Decimal(info["dt_s"])),
        "requested_duration_s": seconds, "dm": info["dm"],
        "sample_sha256": digest.hexdigest(),
        "length_policy": "largest 2/3/5/7-smooth multiple of 16 not exceeding requested length",
    }
    product(destination, ".section.json").write_text(json.dumps(details, indent=2) + "\n")
    print(json.dumps(details, indent=2))


def timed(stage, command, workdir):
    log_path = workdir / (stage + ".log")
    print("Starting {}: {}".format(stage, shlex.join(command) if hasattr(shlex, "join") else " ".join(map(shlex.quote, command))), flush=True)
    started = time.perf_counter()
    with log_path.open("w") as log:
        child = subprocess.Popen(command, cwd=str(workdir), stdout=log, stderr=subprocess.STDOUT)
        _, status, usage = os.wait4(child.pid, 0)
        if os.WIFEXITED(status):
            child.returncode = os.WEXITSTATUS(status)
        else:
            child.returncode = -os.WTERMSIG(status)
    elapsed = time.perf_counter() - started
    result = {
        "stage": stage, "wall_s": elapsed, "user_s": usage.ru_utime,
        "system_s": usage.ru_stime, "cpu_s": usage.ru_utime + usage.ru_stime,
        "max_rss_kib": usage.ru_maxrss / (1024 if sys.platform == "darwin" else 1),
        "returncode": child.returncode, "command": command, "log": str(log_path),
    }
    with (workdir / "metrics.jsonl").open("a") as handle:
        handle.write(json.dumps(result) + "\n")
    print("{}: wall={:.3f}s, CPU={:.3f}s, peak RSS={:.1f} MiB".format(stage, elapsed, result["cpu_s"], result["max_rss_kib"] / 1024), flush=True)
    if child.returncode:
        raise RuntimeError("{} exited with {}; see {}".format(stage, child.returncode, log_path))
    return result


def inventory(data_root):
    trials = {}
    for inf_path in sorted(data_root.glob("DM*_bary/*.inf")):
        dat_path = inf_path.with_suffix(".dat")
        if not dat_path.is_file():
            continue
        info = read_info(inf_path)
        key = int(round(info["dm"] * 100))
        if key not in range(2350, 2551, 10) or abs(info["dm"] - key / 100.0) > 1e-6:
            continue
        if key in trials:
            raise ValueError("More than one input pair at DM {:.2f}".format(key / 100.0))
        if info["barycentric"] != 1 or dat_path.stat().st_size != 4 * info["nsamples"]:
            raise ValueError("Input barycentric flag or byte count is invalid: " + str(inf_path))
        info.update({"dat": str(dat_path.resolve()), "inf": str(inf_path.resolve()),
                     "dat_size_bytes": dat_path.stat().st_size,
                     "dat_mtime_ns": dat_path.stat().st_mtime_ns,
                     "inf_sha256": hashlib.sha256(inf_path.read_bytes()).hexdigest()})
        trials[key] = info
    missing = sorted(set(range(2350, 2551, 10)) - set(trials))
    if missing:
        raise ValueError("Missing DM trials: " + ", ".join("{:.2f}".format(x / 100.0) for x in missing))
    return trials


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Submit run_presto_pilot.sbatch from the host shell to run on a compute node")
    if int(os.environ.get("SLURM_CPUS_PER_TASK", "1")) != 1:
        raise RuntimeError("This pilot expects one allocated CPU per task")
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[name] = "1"
    os.environ["OMP_DYNAMIC"] = "FALSE"
    programs = {}
    for program in ("realfft", "rednoise", "accelsearch"):
        programs[program] = shutil.which(program)
        if not programs[program]:
            raise RuntimeError("Executable is missing from container PATH: " + program)
    trials = inventory(args.data_root.resolve())
    dm_key = int(round(args.dm * 100))
    if dm_key not in trials or abs(args.dm - dm_key / 100.0) > 1e-6:
        raise ValueError("Choose one of the 21 available DM trials")
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    run_dir = args.work_root.resolve() / ("pilot_" + os.environ["SLURM_JOB_ID"] + "_" + stamp)
    if os.path.commonpath([str(run_dir), str(args.data_root.resolve())]) == str(args.data_root.resolve()):
        raise ValueError("Working directory must be outside the source data directory")
    run_dir.mkdir(parents=True, exist_ok=False)
    (run_dir / "input_manifest.json").write_text(json.dumps(list(trials.values()), indent=2) + "\n")
    cpu_lines = []
    if Path("/proc/cpuinfo").exists():
        cpu_lines = [line for line in Path("/proc/cpuinfo").read_text().splitlines() if line.startswith(("model name", "vendor_id", "cpu family", "model\t"))][:4]
    metadata = {
        "status": "running", "run_dir": str(run_dir), "hostname": platform.node(),
        "platform": platform.platform(), "python": sys.version, "cpu": cpu_lines,
        "affinity": sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
        "job_id": os.environ["SLURM_JOB_ID"], "allocated_cpus_per_task": 1,
        "programs": programs, "dm": args.dm, "frequency_hz": [50, 1050],
        "numharm": 1, "zmax": 200, "wmax": 0, "sigma": 2.0,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "image": os.environ.get("PRESTO_PILOT_IMAGE"),
        "candidate_threshold_note": "PRESTO cutoff only; not a matched HMM false-alarm probability",
    }
    summary_path = run_dir / "summary.json"
    summary_path.write_text(json.dumps(metadata, indent=2) + "\n")
    print("Run directory: " + str(run_dir), flush=True)
    source = Path(trials[dm_key]["dat"])
    before_stat = source.stat()
    section = run_dir / "pilot_DM{:.2f}".format(args.dm)
    stages = []
    try:
        stages.append(timed("01_section", [sys.executable, str(Path(__file__).resolve()), "extract", str(source), str(section), "--seconds", str(args.seconds)], run_dir))
        stages.append(timed("02_fft", [programs["realfft"], "-fwd", product(section, ".dat").name], run_dir))
        if not product(section, ".fft").is_file():
            raise RuntimeError("realfft did not produce its expected output")
        stages.append(timed("03_whitening", [programs["rednoise"], product(section, ".fft").name], run_dir))
        red = Path(str(section) + "_red")
        if not product(red, ".fft").is_file() or not product(red, ".inf").is_file():
            raise RuntimeError("rednoise did not produce the expected _red.fft and _red.inf pair")
        stages.append(timed("04_acceleration", [programs["accelsearch"], "-ncpus", "1", "-numharm", "1", "-zmax", "200", "-wmax", "0", "-sigma", "2.0", "-flo", "50", "-fhi", "1050", product(red, ".fft").name], run_dir))
        after_stat = source.stat()
        if (before_stat.st_size, before_stat.st_mtime_ns) != (after_stat.st_size, after_stat.st_mtime_ns):
            raise RuntimeError("Source data size or modification time changed during the pilot")
        metadata.update({"status": "complete", "section": json.loads(product(section, ".section.json").read_text()),
                         "search_wall_s": sum(x["wall_s"] for x in stages[1:]),
                         "search_cpu_s": sum(x["cpu_s"] for x in stages[1:]),
                         "section_plus_search_wall_s": sum(x["wall_s"] for x in stages),
                         "section_plus_search_cpu_s": sum(x["cpu_s"] for x in stages),
                         "peak_tool_rss_kib": max(x["max_rss_kib"] for x in stages),
                         "outputs": sorted(p.name for p in run_dir.iterdir())})
        with (run_dir / "timings.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["stage", "wall_s", "user_s", "system_s", "cpu_s", "max_rss_kib", "returncode"], extrasaction="ignore")
            writer.writeheader()
            writer.writerows(stages)
        print("\n" + (run_dir / "timings.tsv").read_text(), flush=True)
    except Exception as error:
        metadata.update({"status": "failed", "error": str(error)})
        raise
    finally:
        summary_path.write_text(json.dumps(metadata, indent=2) + "\n")
    print("Pilot complete. Summary: " + str(summary_path), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser("run")
    run.add_argument("--data-root", type=Path, default=DATA_ROOT)
    run.add_argument("--work-root", type=Path, default=WORK_ROOT)
    run.add_argument("--dm", type=float, default=24.40)
    run.add_argument("--seconds", type=float, default=900.0)
    cut = commands.add_parser("extract")
    cut.add_argument("source", type=Path)
    cut.add_argument("destination", type=Path)
    cut.add_argument("--seconds", type=float, required=True)
    arguments = parser.parse_args()
    if arguments.command == "extract":
        extract(arguments.source, arguments.destination, arguments.seconds)
    else:
        main(arguments)
