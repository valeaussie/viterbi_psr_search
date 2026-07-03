#!/usr/bin/env python
"""
aggregate_blind.py

Aggregate the blind Viterbi search results across all DM trials, subbands,
and Nt values.  Reads every *_loglike_curve.dat file under the blind output
tree, computes a detection threshold for each distinct Nt value calibrated
to a user-specified false alarm rate per subband, finds peaks above those
thresholds, assembles a master candidate table, and cross-matches candidate
frequencies against the known pulsar catalogue (from Stage 0).

Threshold method (O'Leary, Dunn & Melatos 2026, Section 2.7 and Appendix D)
----------------------------------------------------------------------------
The loglike statistic L = ln P(Q*|O) is a cumulative sum of normalised
periodogram values over N_T coherent segments, so it scales with N_T.
Pooling loglike values across different N_T values gives a multi-modal
distribution that cannot be fitted as a single exponential.  The threshold
is therefore computed separately for each distinct N_T value.

Within each N_T group, all loglike values (one per frequency bin per curve)
from all subbands and DM trials are pooled.  This is equivalent to the
paper's Monte Carlo procedure with N_real = n_curves independent noise
realisations of the spectrogram, where n_curves is the number of curve files
in the group.  Because most subbands are noise-dominated, this approximates
the noise-only distribution.

The procedure follows Eqs 13-17 of the paper:

    lambda_hat = N_tail / sum(L_i - L_tail)                    [Eq. 14]
    L_th = L_tail - (1/lambda_hat)                             [Eq. 17]
              * ln( n_curves * N_Q * (1-(1-alpha')^(1/N_Q)) / N_tail )

where N_Q is the median curve length for the N_T group, L_tail is the
--tail-percentile of the pooled values, and n_curves plays the role of
N_real (number of independent Viterbi runs used to build p(L)).

Deduplication tolerance
-----------------------
The tolerance is derived from the coarsest frequency bin in the Nt grid:

    dedup_tol = k * N_T_max / T_obs

where k = --dedup-bin-widths (default 2) and T_obs is read from the sibling
.params files (T_obs = Tsft * Nt).  Pass --dedup-tol-hz to override.

Two output files are written:

    candidates_raw.csv
    candidates_dedup.csv

Usage
-----
    python aggregate_blind.py \\
        --blind-dir  <stage3_viterbi/blind_v1> \\
        --known-yaml <config/known_pulsars.yaml> \\
        --out-dir    <stage3_viterbi/blind_v1> \\
        [--false-alarm-rate  0.1] \\
        [--tail-percentile   99.99] \\
        [--dedup-bin-widths  2] \\
        [--dedup-tol-hz      <override in Hz>] \\
        [--t-obs             <T_obs in s, fallback if no .params files>] \\
        [--match-tol-hz      0.05]
"""

import argparse
import glob
import os
import re
from collections import defaultdict

import numpy as np

try:
    import yaml
    HAVE_YAML = True
except ImportError:
    HAVE_YAML = False


# ---------------------------------------------------------------------------
# Threshold computation  (paper Eqs 13-17 and Appendix D)
# ---------------------------------------------------------------------------

def compute_paper_threshold(loglike_values, n_curves, n_q, false_alarm_rate,
                             tail_percentile=99.99):
    """
    Compute the detection threshold L_th for one N_T group.

    Parameters
    ----------
    loglike_values : array_like
        All loglike values pooled from curves in this N_T group.
        Total size = n_curves * N_Q (approximately).
    n_curves : int
        Number of independent Viterbi runs (curve files) in this group.
        This is N_real in the paper's notation.
    n_q : int
        Median number of frequency bins per curve in this group (N_Q).
    false_alarm_rate : float
        Desired false alarm probability per subband alpha'.
    tail_percentile : float
        Percentile above which to fit the exponential tail.

    Returns
    -------
    l_tail, lambda_hat, l_th, sigma_l_th, n_tail
    """
    L = np.asarray(loglike_values, dtype=float)
    L = L[np.isfinite(L)]

    if L.size == 0:
        raise ValueError("No finite loglike values.")

    l_tail = float(np.percentile(L, tail_percentile))
    tail = L[L > l_tail]
    n_tail = tail.size

    if n_tail < 10:
        raise ValueError(
            f"Only {n_tail} samples above the {tail_percentile}th percentile."
        )

    # MLE rate for exponential tail  [Eq. 14]
    lambda_hat = float(n_tail / np.sum(tail - l_tail))

    # Per-bin false alarm probability  [Eq. 16 inverted]
    # alpha' = 1 - (1-alpha)^N_Q  =>  alpha = 1 - (1-alpha')^(1/N_Q)
    alpha = 1.0 - (1.0 - false_alarm_rate) ** (1.0 / n_q)

    # L_th  [Eq. 17]
    # n_curves plays the role of N_real: the number of independent
    # Viterbi runs used to construct p(L).
    log_arg = float(n_curves) * float(n_q) * alpha / float(n_tail)
    if log_arg <= 0:
        raise ValueError(f"Invalid log argument {log_arg:.3e}.")
    l_th = l_tail - np.log(log_arg) / lambda_hat

    # Uncertainty via first-order delta method (Appendix D)
    d_lth_d_lambda = np.log(log_arg) / lambda_hat ** 2
    var_lambda = lambda_hat ** 2 / n_tail
    sigma_l_th = float(np.sqrt(d_lth_d_lambda ** 2 * var_lambda))

    return float(l_tail), float(lambda_hat), float(l_th), sigma_l_th, n_tail


# ---------------------------------------------------------------------------
# Deduplication tolerance from bin width
# ---------------------------------------------------------------------------

def read_tsft_from_params(loglike_curve_path):
    """
    Read Tsft from the sibling .params file written by viterbi_pipeline.py.
    Returns float or None.
    """
    params_path = loglike_curve_path.replace("_loglike_curve.dat", ".params")
    if not os.path.isfile(params_path):
        return None
    try:
        with open(params_path) as fh:
            for line in fh:
                line = line.strip()
                if line.startswith("#") or "=" not in line:
                    continue
                key, _, val = line.partition("=")
                if key.strip().lower() == "tsft":
                    return float(val.strip())
    except Exception:
        return None
    return None


def compute_dedup_tolerance(nt_tsft_pairs, n_dedup_bins):
    """
    Derive dedup frequency tolerance as k * N_T_max / T_obs.
    Returns (dedup_tol_hz, t_obs, max_nt).
    """
    t_obs_estimates = [tsft * nt for nt, tsft in nt_tsft_pairs if tsft is not None]
    if not t_obs_estimates:
        raise ValueError(
            "No Tsft values could be read from .params files. "
            "Supply --t-obs or check that .params files exist."
        )
    t_obs = float(np.median(t_obs_estimates))
    max_nt = max(nt for nt, _ in nt_tsft_pairs)
    coarse_bin_hz = max_nt / t_obs
    return n_dedup_bins * coarse_bin_hz, t_obs, max_nt


# ---------------------------------------------------------------------------
# Peak finding
# ---------------------------------------------------------------------------

def find_peaks(freqs, loglike, loglike_threshold, min_separation_bins=5):
    """
    Return list of (peak_freq_hz, peak_loglike, threshold) tuples for all
    local maxima that exceed loglike_threshold.
    """
    L = np.asarray(loglike, dtype=float)
    f = np.asarray(freqs, dtype=float)
    if L.size < 3:
        return []

    peaks = []
    last_peak_idx = -(10 ** 9)
    for i in range(1, L.size - 1):
        if L[i] >= L[i - 1] and L[i] > L[i + 1] and L[i] > loglike_threshold:
            if i - last_peak_idx >= min_separation_bins:
                peaks.append((float(f[i]), float(L[i]), float(loglike_threshold)))
                last_peak_idx = i
            elif peaks and L[i] > peaks[-1][1]:
                peaks[-1] = (float(f[i]), float(L[i]), float(loglike_threshold))
                last_peak_idx = i
    return peaks


# ---------------------------------------------------------------------------
# Path metadata parsing
# ---------------------------------------------------------------------------

def parse_run_metadata(path):
    """
    Extract DM, Nt, and subband f0 from a path of the form:
        .../blind_v1/DM<XX.XX>/Nt<N>/f0_<F>/blind_Nt<N>_f0_<F>_loglike_curve.dat
    Returns (dm, nt, f0_sub). Unparseable fields returned as None.
    """
    dm = nt = f0_sub = None

    m = re.search(r"/DM([0-9]+\.[0-9]+)/", path)
    if m:
        dm = float(m.group(1))

    m = re.search(r"/Nt(\d+)/", path)
    if m:
        nt = int(m.group(1))

    m = re.search(r"/f0_([0-9.]+)/", path)
    if m:
        f0_sub = float(m.group(1))

    return dm, nt, f0_sub


# ---------------------------------------------------------------------------
# Known pulsar catalogue
# ---------------------------------------------------------------------------

def load_known_pulsars(yaml_path):
    """Return list of (name, F0_hz) for known pulsars with an F0."""
    if not (HAVE_YAML and yaml_path and os.path.exists(yaml_path)):
        return []
    with open(yaml_path) as fh:
        data = yaml.safe_load(fh)
    out = []
    for p in data.get("pulsars", []):
        f0 = p.get("F0_hz")
        name = p.get("name", "?")
        if f0 is not None:
            out.append((name, float(f0)))
    return out


def match_known(freq, known, tol_hz):
    """Return name of known pulsar within tol_hz of freq, else empty string."""
    best = ""
    best_df = tol_hz
    for name, f0 in known:
        df = abs(freq - f0)
        if df <= best_df:
            best = name
            best_df = df
    return best


# ---------------------------------------------------------------------------
# Deduplication
# ---------------------------------------------------------------------------

def dedup_candidates(rows, freq_tol_hz):
    """
    Collapse peaks within freq_tol_hz of each other into single candidates,
    merging across DM trials, Nt values, and overlapping subbands.
    Returns list of dicts sorted by peak_loglike descending.
    """
    ordered = sorted(rows, key=lambda r: r["peak_loglike"], reverse=True)
    clusters = []

    for r in ordered:
        f = r["peak_freq_hz"]
        placed = False
        for c in clusters:
            if abs(f - c["peak_freq_hz"]) <= freq_tol_hz:
                c["multiplicity"] += 1
                c["nt_set"].add(r["Nt"])
                c["dm_set"].add(r["dm"])
                if not c["known_match"] and r["known_match"]:
                    c["known_match"] = r["known_match"]
                placed = True
                break
        if not placed:
            clusters.append({
                "peak_freq_hz": f,
                "peak_loglike": r["peak_loglike"],
                "dm_best":      r["dm"],
                "Nt_best":      r["Nt"],
                "subband_f0":   r["subband_f0"],
                "threshold":    r["threshold"],
                "known_match":  r["known_match"],
                "multiplicity": 1,
                "nt_set":       {r["Nt"]},
                "dm_set":       {r["dm"]},
            })

    for c in clusters:
        c["nt_values"] = ";".join(
            str(x) for x in sorted(v for v in c["nt_set"] if v is not None))
        c["dm_values"] = ";".join(
            f"{v:.2f}" for v in sorted(v for v in c["dm_set"] if v is not None))
        c["dm_count"] = len(c["dm_set"])
        del c["nt_set"]
        del c["dm_set"]

    clusters.sort(key=lambda c: c["peak_loglike"], reverse=True)
    return clusters


# ---------------------------------------------------------------------------
# CSV helpers
# ---------------------------------------------------------------------------

def write_csv(path, cols, rows):
    with open(path, "w") as fh:
        fh.write(",".join(cols) + "\n")
        for r in rows:
            fh.write(",".join(str(r.get(c, "")) for c in cols) + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Aggregate blind Viterbi search results across DM, Nt, and "
                    "subband grid using the paper's exponential-tail threshold.")
    ap.add_argument("--blind-dir", required=True)
    ap.add_argument("--known-yaml", default=None)
    ap.add_argument("--out-dir", required=True)

    # Threshold
    ap.add_argument("--false-alarm-rate", type=float, default=0.1,
                    help="False alarm probability per subband alpha' (default 0.1).")
    ap.add_argument("--tail-percentile", type=float, default=99.99,
                    help="Percentile for exponential tail cut-off (default 99.99).")

    # Dedup
    ap.add_argument("--dedup-bin-widths", type=float, default=2.0,
                    help="Dedup tolerance as multiple of coarsest bin width (default 2).")
    ap.add_argument("--dedup-tol-hz", type=float, default=None,
                    help="Override: fixed dedup tolerance in Hz.")
    ap.add_argument("--t-obs", type=float, default=None,
                    help="T_obs in seconds, fallback if .params files are unreadable.")

    # Cross-match
    ap.add_argument("--match-tol-hz", type=float, default=0.05)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    pattern = os.path.join(args.blind_dir, "**", "*_loglike_curve.dat")
    files = sorted(glob.glob(pattern, recursive=True))
    if not files:
        print(f"No *_loglike_curve.dat files found under {args.blind_dir}")
        return

    known = load_known_pulsars(args.known_yaml)
    print(f"Loaded {len(known)} known pulsars with F0.")
    print(f"Found {len(files)} loglike curve files.")

    # ------------------------------------------------------------------
    # Pass 1: load all curves, group by Nt
    # ------------------------------------------------------------------
    print("\nPass 1: loading curves and grouping by Nt...")

    groups = defaultdict(list)   # nt -> list of (path, freqs, loglike, dm, f0_sub)
    nt_tsft_pairs = []
    n_skipped = 0

    for path in files:
        try:
            arr = np.loadtxt(path)
        except Exception as exc:
            print(f"  skip {path}: {exc}")
            n_skipped += 1
            continue
        if arr.ndim != 2 or arr.shape[0] < 3:
            n_skipped += 1
            continue

        freqs, loglike = arr[:, 0], arr[:, 1]
        dm, nt, f0_sub = parse_run_metadata(path)
        tsft = read_tsft_from_params(path)

        if nt is None:
            n_skipped += 1
            continue

        groups[nt].append((path, freqs, loglike, dm, f0_sub))
        nt_tsft_pairs.append((nt, tsft))

    if not groups:
        print("ERROR: no usable loglike curves found.")
        return

    if n_skipped:
        print(f"  ({n_skipped} files skipped)")

    nt_values_found = sorted(groups.keys())
    total_curves = sum(len(v) for v in groups.values())
    print(f"  {total_curves} curves across Nt = {nt_values_found}")

    # ------------------------------------------------------------------
    # Deduplication tolerance
    # ------------------------------------------------------------------
    if args.dedup_tol_hz is not None:
        dedup_tol_hz = args.dedup_tol_hz
        dedup_source = "user override (--dedup-tol-hz)"
        t_obs = None
        max_nt = None
    else:
        if args.t_obs is not None:
            nt_tsft_pairs = [
                (nt, tsft if tsft is not None else args.t_obs / nt)
                for nt, tsft in nt_tsft_pairs
            ]
        try:
            dedup_tol_hz, t_obs, max_nt = compute_dedup_tolerance(
                nt_tsft_pairs, args.dedup_bin_widths
            )
            coarse_bin_hz = max_nt / t_obs
            dedup_source = (
                f"auto: {args.dedup_bin_widths} x coarsest bin "
                f"({args.dedup_bin_widths} x {coarse_bin_hz:.5f} Hz, "
                f"N_T_max={max_nt}, T_obs={t_obs:.1f} s)"
            )
        except ValueError as exc:
            print(f"\nERROR computing dedup tolerance: {exc}")
            return

    print(f"\nDeduplication tolerance: {dedup_tol_hz:.5f} Hz  [{dedup_source}]")

    # ------------------------------------------------------------------
    # Pass 2: per-Nt threshold calibration and peak finding
    # ------------------------------------------------------------------
    print("\nPass 2: per-Nt threshold calibration and peak finding...")
    print(f"  (N_real for each Nt group = number of curves in that group)")

    all_rows = []
    threshold_records = []

    for nt in nt_values_found:
        curves = groups[nt]
        n_curves_nt = len(curves)

        all_loglike = np.concatenate([loglike for _, _, loglike, _, _ in curves])
        n_q_median = int(np.median([len(loglike) for _, _, loglike, _, _ in curves]))

        try:
            l_tail, lambda_hat, l_th, sigma_l_th, n_tail = compute_paper_threshold(
                all_loglike,
                n_curves=n_curves_nt,
                n_q=n_q_median,
                false_alarm_rate=args.false_alarm_rate,
                tail_percentile=args.tail_percentile,
            )
        except ValueError as exc:
            print(f"  Nt={nt}: skipping threshold ({exc})")
            continue

        print(f"  Nt={nt:>4d}: {n_curves_nt:>5d} curves  "
              f"N_Q={n_q_median}  "
              f"L_tail={l_tail:.2f}  "
              f"lambda={lambda_hat:.4f}  "
              f"L_th={l_th:.2f} +/- {sigma_l_th:.2f}")

        threshold_records.append({
            "nt": nt, "n_curves": n_curves_nt, "n_q_median": n_q_median,
            "l_tail": l_tail, "n_tail": n_tail, "lambda_hat": lambda_hat,
            "l_th": l_th, "sigma_l_th": sigma_l_th,
        })

        for path, freqs, loglike, dm, f0_sub in curves:
            peaks = find_peaks(freqs, loglike, l_th)
            for (pf, pl, thr) in peaks:
                match = match_known(pf, known, args.match_tol_hz)
                all_rows.append({
                    "dm":           dm,
                    "Nt":           nt,
                    "subband_f0":   f0_sub,
                    "peak_freq_hz": pf,
                    "peak_loglike": pl,
                    "threshold":    thr,
                    "known_match":  match,
                })

    all_rows.sort(key=lambda r: r["peak_loglike"], reverse=True)

    # ------------------------------------------------------------------
    # Write candidates_raw.csv
    # ------------------------------------------------------------------
    raw_cols = ["dm", "Nt", "subband_f0", "peak_freq_hz",
                "peak_loglike", "threshold", "known_match"]
    raw_csv = os.path.join(args.out_dir, "candidates_raw.csv")
    write_csv(raw_csv, raw_cols, all_rows)

    n_known_raw = sum(1 for r in all_rows if r["known_match"])
    print(f"\nWrote {raw_csv}")
    print(f"  total peaks above threshold : {len(all_rows)}")
    print(f"  matched to known pulsars    : {n_known_raw}")
    print(f"  unmatched                   : {len(all_rows) - n_known_raw}")

    # ------------------------------------------------------------------
    # Write candidates_dedup.csv
    # ------------------------------------------------------------------
    clusters = dedup_candidates(all_rows, dedup_tol_hz)

    dedup_cols = [
        "peak_freq_hz", "peak_loglike",
        "dm_best", "dm_count", "dm_values",
        "multiplicity", "nt_values", "Nt_best",
        "subband_f0", "threshold", "known_match",
    ]
    dedup_csv = os.path.join(args.out_dir, "candidates_dedup.csv")
    write_csv(dedup_csv, dedup_cols, clusters)

    n_known_d = sum(1 for c in clusters if c["known_match"])
    known_names = sorted({c["known_match"] for c in clusters if c["known_match"]})
    print(f"\nWrote {dedup_csv}")
    print(f"  distinct candidates (tol={dedup_tol_hz:.5f} Hz) : {len(clusters)}")
    print(f"  distinct known pulsars recovered               : {len(known_names)}")
    if known_names:
        print(f"    {', '.join(known_names)}")
    print(f"  distinct unmatched candidates                  : "
          f"{len(clusters) - n_known_d}")

    # ------------------------------------------------------------------
    # Write calibration metadata
    # ------------------------------------------------------------------
    threshold_path = os.path.join(args.out_dir, "threshold_calibration.txt")
    with open(threshold_path, "w") as fh:
        fh.write("# Threshold calibration per Nt  (O'Leary, Dunn & Melatos 2026)\n")
        fh.write(f"false_alarm_rate  = {args.false_alarm_rate}\n")
        fh.write(f"tail_percentile   = {args.tail_percentile}\n")
        fh.write(f"dedup_tol_hz      = {dedup_tol_hz:.6f}\n")
        fh.write(f"dedup_source      = {dedup_source}\n")
        if t_obs is not None:
            fh.write(f"t_obs_inferred_s  = {t_obs:.3f}\n")
        if max_nt is not None:
            fh.write(f"max_nt            = {max_nt}\n")
        fh.write("\n")
        fh.write(f"{'nt':>6}  {'n_curves':>8}  {'n_q':>6}  {'l_tail':>10}  "
                 f"{'n_tail':>8}  {'lambda':>8}  {'l_th':>10}  {'sigma':>8}\n")
        for r in threshold_records:
            fh.write(
                f"{r['nt']:>6d}  {r['n_curves']:>8d}  {r['n_q_median']:>6d}  "
                f"{r['l_tail']:>10.4f}  {r['n_tail']:>8d}  "
                f"{r['lambda_hat']:>8.4f}  {r['l_th']:>10.4f}  "
                f"{r['sigma_l_th']:>8.4f}\n"
            )
    print(f"\nCalibration metadata written to {threshold_path}")

    print("\nTop 20 distinct candidates:")
    hdr = (f"  {'freq_hz':>12}  {'loglike':>9}  {'dm_best':>7}  "
           f"{'dm_count':>8}  {'mult':>5}  {'nt_values':>16}  known_match")
    print(hdr)
    for c in clusters[:20]:
        print(f"  {c['peak_freq_hz']:>12.4f}  {c['peak_loglike']:>9.1f}  "
              f"{c['dm_best']:>7.2f}  {c['dm_count']:>8}  "
              f"{c['multiplicity']:>5}  {c['nt_values']:>16}  {c['known_match']}")


if __name__ == "__main__":
    main()