"""Plot all selected HMM candidates using saved arrays and saved fits only."""
from pathlib import Path
from decimal import Decimal, getcontext
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import tarfile
import time

ROOT = Path(__file__).resolve().parent
BASE = ROOT.parent
INPUT = ROOT / "inputs"
CONFIG = ROOT / "qa" / "matplotlib"
os.environ["MPLCONFIGDIR"] = str(CONFIG)
os.environ["MPLBACKEND"] = "Agg"
getcontext().prec = 30

def extract_inputs():
    marker = INPUT / "extracted.json"
    if marker.exists():
        return
    with tarfile.open(INPUT / "saved_inputs.tar") as archive:
        members = archive.getmembers()
        for member in members:
            target = INPUT / member.name
            if not member.isfile() or not target.resolve().is_relative_to(INPUT.resolve()):
                raise ValueError("Unexpected archive entry: " + member.name)
        archive.extractall(INPUT, filter="data")
    marker.write_text(json.dumps(dict(files=len(members))) + "\n")

def summary(meta, models, epoch_string, df, nt):
    poly, kep = models["polynomial"], models["circular"]
    epoch = Decimal(epoch_string)
    poly_epoch = epoch + Decimal(str(poly["reference_time_s"])) / Decimal(86400)
    t0_epoch = epoch + Decimal(str(kep["parameters"][3])) / Decimal(86400)
    lines = [
        "# Viterbi-search fit summary: original selected track and saved fits",
        "# No search, refitting, frequency refinement or folding was performed.",
        "candidate_id " + meta["candidate_id"],
        "classification " + (meta.get("classification") or "unclassified"),
        "parent_pulsar " + (meta.get("parent_pulsar") or "none"),
        "rfi_flag " + (meta.get("rfi_flag") or "none"),
        "selected_candidate_retained True",
        "frequency_frame barycentric",
        "time_reference seconds_from_barycentric_input_epoch",
        "input_epoch_mjd " + epoch_string,
        "# The .inf file identifies barycentering but does not state the time scale.",
        f"dm_pc_cm-3 {float(meta['dm_best']):.2f}",
        f"n_segments {nt}",
        f"segment_duration_s {meta['sampling']['segment_duration_s']:.15g}",
        f"search_frequency_bin_hz {df:.15e}",
        f"selected_ending_frequency_hz {float(meta['peak_freq_hz']):.15f}",
        f"selected_path_score {float(meta['peak_loglike']):.15f}",
        f"selection_threshold {float(meta['threshold']):.15f}",
        "# The blind threshold is a ranking prescription, not a calibrated false-alarm probability.",
        "", "# Polynomial fit: f(t) = sum_k f_k (t-tref)^k / k!",
        "# Columns: parameter, saved fitted value, saved formal 1-sigma error.",
        "# f_k units are Hz s^(-k); f_0 is the apparent frequency at tref.",
        f"poly_order_requested {poly['order']}",
        f"poly_order_used {poly['order']}",
        f"poly_tref_s {poly['reference_time_s']:.12f}",
        f"poly_tref_mjd {poly_epoch:.15f}",
        f"# tref is the centre of segment index {nt // 2} (zero-based), not the observation midpoint.",
    ]
    for k, (v, e) in enumerate(zip(poly["frequency_derivatives"], poly["frequency_derivative_errors"])):
        lines.append(f"poly_f{k} {v:.15e} {e:.15e}")
    lines += [f"poly_chi2 {poly['chi2']:.15g}", f"poly_dof {poly['dof']}",
              f"poly_red_chi2 {poly['reduced_chi2']:.15g}", "",
              "# Circular Kepler fit: f_app(t) = f_spin - K sin(2*pi*(t-T0)/Pb)",
              "# f_spin is the fitted model's constant frequency; K is in Hz.",
              "# Pb and T0 are in seconds; T0 is relative to the barycentric input epoch.",
              "# T0 is the model's downward zero crossing, not a verified T_asc convention."]
    names = ("kepler_f_spin_hz", "kepler_K_hz", "kepler_Pb_s", "kepler_T0_s")
    for name, value, error in zip(names, kep["parameters"], kep["parameter_errors"]):
        lines.append(f"{name} {value:.15e} {error:.15e}")
    lines += [f"kepler_T0_mjd {t0_epoch:.15f}", f"kepler_chi2 {kep['chi2']:.15g}",
              f"kepler_dof {kep['dof']}", f"kepler_red_chi2 {kep['reduced_chi2']:.15g}", "",
              "# Both fits describe a quantized frequency track, not a phase-connected timing solution.",
              "# Saved errors are fit-covariance estimates, not calibrated detection significances.",
              "# No pulsar identification or orbital detection is inferred from these fits.",
              "# The complete circular model and all six polynomial derivatives are retained."]
    return "\n".join(lines) + "\n"

def process_source(task):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator, ScalarFormatter

    source_index, metas = task
    cached = INPUT / "source_cache" / f"source_{source_index:03d}"
    array_path = cached / "original_search_state.npz"
    curve_path = cached / "selected_source_curve.dat"
    with np.load(array_path, allow_pickle=False) as state:
        freq = state["frequency_hz"]
        spec = state["spectrogram"]
        delta = state["delta"]
        tsft = float(state["segment_duration_s"])
    curve = np.loadtxt(curve_path)
    nf, nt = spec.shape
    df = float(np.median(np.diff(freq)))
    duration = nt * tsft
    assert delta.shape == spec.shape and nf == len(freq)
    np.testing.assert_array_equal(curve[:, 0], freq)
    np.testing.assert_allclose(curve[:, 1], delta[:, -1], rtol=0, atol=1e-6)
    source_hash = hashlib.sha256(array_path.read_bytes()).hexdigest()
    curve_hash = hashlib.sha256(curve_path.read_bytes()).hexdigest()
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.labelsize": 11, "axes.titlesize": 10, "legend.fontsize": 8.5,
        "pdf.fonttype": 42, "axes.linewidth": 0.8,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
    })
    reports = []
    script_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    for meta in metas:
        cid = meta["candidate_id"]
        output = ROOT / "candidates" / cid
        complete = output / "verification.json"
        if complete.exists():
            report = json.loads(complete.read_text())
            if report["plotting_script_sha256"] == script_hash:
                reports.append(report)
                continue
        original = BASE / "hmm_gallery" / cid
        saved = BASE / "fold_gallery" / "hmm" / cid
        models = json.loads((saved / "models.json").read_text())
        fitdata = np.loadtxt(saved / "track_with_fits.dat")
        track = np.loadtxt(original / (cid + "_track.dat"))
        t, path = track.T
        poly, kep = models["polynomial"], models["circular"]
        fspin, amplitude, pb, t0 = kep["parameters"]
        end = int(meta["source_bin"])
        threshold = float(meta["threshold"])
        peak_frequency = float(meta["peak_freq_hz"])
        peak_score = float(meta["peak_loglike"])
        assert nt == meta["sampling"]["actual_segments"]
        assert meta["exact_original_path"] and meta["selected_candidate_retained"]
        assert tsft == meta["sampling"]["segment_duration_s"]
        np.testing.assert_array_equal(track, fitdata[:, :2])
        np.testing.assert_allclose(t, (np.arange(nt) + .5) * tsft, rtol=0, atol=1e-9)
        indices = np.empty(nt, dtype=int)
        indices[-1] = end
        for j in range(nt - 1, 0, -1):
            i = indices[j]
            a, b = max(0, i - 1), min(nf, i + 2)
            indices[j - 1] = a + np.argmax(delta[a:b, j - 1])
        np.testing.assert_array_equal(path, freq[indices])
        assert abs(freq[end] - peak_frequency) < 1e-9
        assert abs(curve[end, 1] - peak_score) < 1e-6

        def polynomial(times):
            return np.polynomial.polynomial.polyval(
                np.asarray(times) - poly["reference_time_s"],
                [v / math.factorial(k) for k, v in enumerate(poly["frequency_derivatives"])])

        def circular(times):
            return fspin - amplitude * np.sin(2 * np.pi * (np.asarray(times) - t0) / pb)

        np.testing.assert_allclose(polynomial(t), fitdata[:, 2], rtol=0, atol=1e-9)
        np.testing.assert_allclose(circular(t), fitdata[:, 3], rtol=0, atol=1e-9)
        epoch_path = INPUT / "epochs" / f"47Tuc_DM{float(meta['dm_best']):.2f}_bary.inf"
        inf = {}
        for line in epoch_path.read_text().splitlines():
            if "=" in line:
                k, v = line.split("=", 1)
                inf[k.strip()] = v.strip()
        assert next(v for k, v in inf.items() if k.startswith("Barycentered?")) == "1"
        assert float(inf["Dispersion measure (cm-3 pc)"]) == float(meta["dm_best"])
        epoch_string = inf["Epoch of observation (MJD)"]
        output.mkdir(parents=True, exist_ok=True)
        prefix = output / cid

        def save(fig, suffix):
            fig.savefig(str(prefix) + suffix + ".pdf", dpi=300)
            fig.savefig(str(prefix) + suffix + ".png", dpi=160)
            plt.close(fig)

        title = (f"Candidate {int(meta['cand_index']):06d}: "
                 + rf"DM = {float(meta['dm_best']):.2f} pc cm$^{{-3}}$, $N_T = {nt}$")
        identity = meta.get("parent_pulsar") or "unclassified"
        flag = meta.get("rfi_flag") or "none"
        subtitle = f"Barycentric; {identity}; RFI flag: {flag} (retained)"
        dense = np.unique(np.r_[np.linspace(0, duration, max(2001, int(64 * duration / pb))), t])
        poly_dense, kep_dense = polynomial(dense), circular(dense)
        low = min(path.min(), poly_dense.min(), kep_dense.min())
        high = max(path.max(), poly_dense.max(), kep_dense.max())
        margin = max(5 * df, .06 * (high - low))
        lo, hi = max(freq[0] - df / 2, low - margin), min(freq[-1] + df / 2, high + margin)
        mask = (freq >= lo) & (freq <= hi)
        view = spec[mask]
        fig, ax = plt.subplots(figsize=(7.4, 4.9), layout="constrained")
        ax.imshow(view, origin="lower", aspect="auto", interpolation="nearest",
                  extent=(0, duration, freq[mask][0] - df / 2, freq[mask][-1] + df / 2),
                  cmap="gray_r", vmin=np.percentile(view, 20), vmax=np.percentile(view, 99.5), rasterized=True)
        ax.stairs(path, np.arange(nt + 1) * tsft, baseline=None, color="red", lw=1.6,
                  label=rf"Viterbi path ($T_{{\rm coh}} = {tsft:.2f}$ s)")
        ax.plot(dense, poly_dense, color="tab:blue", ls="--", lw=1.6,
                label=f"Polynomial fit (order {poly['order']})")
        ax.plot(dense, kep_dense, color="tab:orange", ls=":", lw=1.9, label="Circular Kepler fit")
        ax.set(xlim=(0, duration), ylim=(lo, hi), xlabel="Time from barycentric input epoch (s)",
               ylabel="Barycentric frequency (Hz)", title=title + "\n" + subtitle)
        ax.ticklabel_format(axis="y", style="plain", useOffset=False)
        ax.legend(loc="upper right", framealpha=.95)
        save(fig, "_spectrogram")

        fig, ax = plt.subplots(figsize=(7.4, 4.7), layout="constrained")
        ax.plot(curve[:, 0], curve[:, 1], color="black", lw=.65, label="Best-path ending-state score")
        ax.axhline(threshold, color="tab:blue", ls="--", lw=1.1,
                   label=f"Selection threshold = {threshold:.4f}")
        ax.plot(freq[end], curve[end, 1], "o", color="red", ms=4, label=f"Selected endpoint = {peak_score:.4f}")
        ax.set(xlim=(freq[0], freq[-1]), xlabel="Final-state barycentric frequency (Hz)",
               ylabel="Accumulated path score", title=title + "\n" + subtitle)
        ax.ticklabel_format(axis="x", style="plain", useOffset=False)
        # Keep the selected peak clear; put the inset and legend on the other side.
        left = (peak_frequency - freq[0]) / (freq[-1] - freq[0]) >= .5
        ax.legend(loc="upper left" if left else "upper right", framealpha=.95)
        inset = ax.inset_axes([.05 if left else .57, .38, .38, .40])
        halfwidth = max(.03, 10 * df)
        pmask = np.abs(freq - peak_frequency) <= halfwidth
        inset.plot(curve[pmask, 0], curve[pmask, 1], color="black", lw=.8)
        inset.axhline(threshold, color="tab:blue", ls="--", lw=.9)
        inset.plot(freq[end], curve[end, 1], "o", color="red", ms=3)
        inset.set_title("Selected endpoint", fontsize=8)
        inset.tick_params(labelsize=7)
        formatter = ScalarFormatter(useOffset=False)
        formatter.set_scientific(False)
        inset.xaxis.set_major_formatter(formatter)
        inset.xaxis.set_major_locator(MaxNLocator(nbins=3))
        save(fig, "_loglikes")

        # User's additional plain blue curve: 2 Hz around the selected endpoint,
        # using the ORIGINAL scores only. Clip the display at the source-band edges.
        window = (freq >= peak_frequency - 1.) & (freq <= peak_frequency + 1.)
        fig, ax = plt.subplots(figsize=(7.4, 4.7), layout="constrained")
        ax.plot(curve[window, 0], curve[window, 1], color="tab:blue", lw=1.5)
        ax.set(xlabel="Frequency (Hz)", ylabel="Log likelihood", title=title + "\n" + subtitle)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.ticklabel_format(axis="x", style="plain", useOffset=False)
        save(fig, "_loglikes_plain_2Hz")

        Path(str(prefix) + "_fit_summary.dat").write_text(summary(meta, models, epoch_string, df, nt))
        shutil.copyfile(saved / "models.json", output / "models.json")
        shutil.copyfile(saved / "track_with_fits.dat", output / "track_with_fits.dat")
        # Full-band numerical likelihoods keep every original ending-state score.
        shutil.copyfile(curve_path, Path(str(prefix) + "_loglikes.dat"))
        report = dict(candidate_id=cid, cand_index=int(meta["cand_index"]), dm=float(meta["dm_best"]),
                      nt=nt, classification=meta.get("classification", "unclassified"),
                      parent_pulsar=meta.get("parent_pulsar", ""), rfi_flag=meta.get("rfi_flag", ""),
                      retained=True, endpoint_hz=peak_frequency, score=peak_score, threshold=threshold,
                      source_index=source_index, original_track_verified=True, full_score_curve_verified=True,
                      saved_fits_verified=True, input_epoch_mjd=epoch_string,
                      polynomial_max_difference_hz=float(np.max(np.abs(polynomial(t) - fitdata[:, 2]))),
                      circular_max_difference_hz=float(np.max(np.abs(circular(t) - fitdata[:, 3]))),
                      model_outside_search_band=bool(low < freq[0] - df / 2 or high > freq[-1] + df / 2),
                      spectrogram_frequency_limits_hz=[lo, hi],
                      plain_likelihood_frequency_limits_hz=[float(freq[window][0]), float(freq[window][-1])],
                      array_sha256=source_hash, curve_sha256=curve_hash,
                      model_sha256=hashlib.sha256((saved / "models.json").read_bytes()).hexdigest(),
                      epoch_sha256=hashlib.sha256(epoch_path.read_bytes()).hexdigest(),
                      plotting_script_sha256=script_hash, search_rerun=False, refitted=False)
        complete.write_text(json.dumps(report, indent=2) + "\n")
        reports.append(report)
    return reports

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--candidates", help="Optional comma-separated candidate indices for preview")
    args = parser.parse_args()
    extract_inputs()
    CONFIG.mkdir(parents=True, exist_ok=True)
    fontcache = BASE / "restored_example_000000/qa/matplotlib/fontlist-v390.json"
    if fontcache.is_file() and not (CONFIG / fontcache.name).exists():
        shutil.copyfile(fontcache, CONFIG / fontcache.name)
    records = sorted((json.loads(p.read_text()) for p in (BASE / "hmm_gallery").glob("*/candidate.json")),
                     key=lambda r: int(r["cand_index"]))
    assert len(records) == 415 and len({r["candidate_id"] for r in records}) == 415
    selected = set(map(int, args.candidates.split(","))) if args.candidates else None
    groups = {}
    for record in records:
        if selected is None or int(record["cand_index"]) in selected:
            groups.setdefault(int(record["source_recovery_index"]), []).append(record)
    started = time.perf_counter()
    finished = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(process_source, item): item[0] for item in sorted(groups.items())}
        for future in as_completed(futures):
            reports = future.result()
            finished.extend(reports)
            print(json.dumps(dict(done=len(finished), expected=sum(map(len, groups.values())),
                                  source=futures[future], elapsed_s=round(time.perf_counter()-started, 1))), flush=True)
    print("Finished approved saved-data plotting", len(finished), flush=True)

if __name__ == "__main__":
    main()
