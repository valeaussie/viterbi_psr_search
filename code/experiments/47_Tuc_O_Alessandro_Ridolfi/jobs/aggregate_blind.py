#!/usr/bin/env python
"""
aggregate_blind.py

Aggregate the blind Viterbi search results across all DM trials, subbands,
and Nt values.  Reads every *_loglike_curve.dat file under the blind output
tree, finds peaks in each loglike-vs-frequency curve, assembles a master
candidate table, and cross-matches candidate frequencies against the known
pulsar catalogue (from Stage 0).

Peak threshold (per curve):
    L_peak > median(L) + n_sigma * 1.4826 * MAD(L)

Two output files are written:

    candidates_raw.csv
        One row per (DM, Nt, subband, peak_frequency) detection.
        Full information, nothing collapsed.

    candidates_dedup.csv
        Deduplicated in frequency only (within --dedup-tol-hz, across
        subbands and Nt values).  DM information is preserved as new
        columns: dm_best, dm_values, dm_count.

Usage
-----
    python aggregate_blind.py \\
        --blind-dir  <stage3_viterbi/blind_v1> \\
        --known-yaml <config/known_47tuc_pulsars.yaml> \\
        --out-dir    <stage3_viterbi/blind_v1> \\
        [--n-sigma 8] [--min-loglike 0] \\
        [--match-tol-hz 0.05] [--dedup-tol-hz 0.5]
"""

import argparse
import glob
import os
import re

import numpy as np

try:
    import yaml
    HAVE_YAML = True
except ImportError:
    HAVE_YAML = False


# ---------------------------------------------------------------------------
# Peak finding
# ---------------------------------------------------------------------------

def find_peaks(freqs, loglike, n_sigma, min_loglike, min_separation_bins=5):
    """Return list of (peak_freq, peak_loglike, threshold) tuples."""
    L = np.asarray(loglike, dtype=float)
    f = np.asarray(freqs, dtype=float)
    if L.size < 3:
        return []

    med = np.median(L)
    mad = np.median(np.abs(L - med))
    sigma = 1.4826 * mad if mad > 0 else np.std(L)
    thr = max(med + n_sigma * sigma, min_loglike)

    peaks = []
    last_peak_idx = -(10 ** 9)
    for i in range(1, L.size - 1):
        if L[i] >= L[i - 1] and L[i] > L[i + 1] and L[i] > thr:
            if i - last_peak_idx >= min_separation_bins:
                peaks.append((float(f[i]), float(L[i]), float(thr)))
                last_peak_idx = i
            elif peaks and L[i] > peaks[-1][1]:
                peaks[-1] = (float(f[i]), float(L[i]), float(thr))
                last_peak_idx = i
    return peaks


# ---------------------------------------------------------------------------
# Path metadata parsing
# ---------------------------------------------------------------------------

def parse_run_metadata(path):
    """
    Extract DM, Nt, and subband f0 from a path of the form:
        .../blind_v1/DM<XX.XX>/Nt<N>/f0_<F>/blind_Nt<N>_f0_<F>_loglike_curve.dat
    Returns (dm, nt, f0_sub) with dm as float, nt as int, f0_sub as float.
    Any field that cannot be parsed is returned as None.
    """
    dm = None
    nt = None
    f0_sub = None

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
# Deduplication (frequency only; DM info preserved as columns)
# ---------------------------------------------------------------------------

def dedup_candidates(rows, freq_tol_hz):
    """
    Collapse peaks within freq_tol_hz of each other into single distinct
    candidates, merging across DM trials, Nt values, and overlapping subbands.

    Greedy: process peaks in descending loglike order.  Each peak either
    joins an existing cluster (if within freq_tol_hz of its representative
    frequency) or starts a new one.

    Each cluster records:
        peak_freq_hz    representative frequency (from highest-loglike peak)
        peak_loglike    maximum loglike across all merges
        dm_best         DM of the highest-loglike detection
        dm_values       sorted CSV of distinct DMs that recovered the candidate
        dm_count        number of distinct DM trials that recovered it
        Nt_best         Nt of the highest-loglike detection
        subband_f0      subband lower edge of the highest-loglike detection
        multiplicity    total number of raw peaks merged (across DM, Nt, subband)
        nt_values       sorted CSV of distinct Nt values that recovered it
        threshold       threshold of the highest-loglike detection
        known_match     matched known pulsar name if any

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
        description="Aggregate blind Viterbi search results across DM, Nt, and subband grid.")
    ap.add_argument("--blind-dir", required=True,
                    help="Root of the blind output tree "
                         "(contains DM*/Nt*/f0_*/ subdirectories).")
    ap.add_argument("--known-yaml", default=None,
                    help="Path to known_pulsars.yaml for cross-matching.")
    ap.add_argument("--out-dir", required=True,
                    help="Directory for output CSV files.")
    ap.add_argument("--n-sigma", type=float, default=8.0,
                    help="Robust peak threshold in MAD-sigma above median (default 8).")
    ap.add_argument("--min-loglike", type=float, default=0.0,
                    help="Absolute minimum loglike for a peak (default 0).")
    ap.add_argument("--match-tol-hz", type=float, default=0.05,
                    help="Frequency tolerance for known-pulsar cross-match (default 0.05 Hz).")
    ap.add_argument("--dedup-tol-hz", type=float, default=0.5,
                    help="Frequency tolerance for deduplication across Nt, subband, "
                         "and DM (default 0.5 Hz).")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    pattern = os.path.join(args.blind_dir, "**", "*_loglike_curve.dat")
    files = sorted(glob.glob(pattern, recursive=True))
    if not files:
        print(f"No *_loglike_curve.dat files found under {args.blind_dir}")
        return

    known = load_known_pulsars(args.known_yaml)
    print(f"Loaded {len(known)} known pulsars with F0.")
    print(f"Scanning {len(files)} loglike curves...")

    rows = []
    n_skipped = 0
    for path in files:
        dm, nt, f0_sub = parse_run_metadata(path)
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
        peaks = find_peaks(freqs, loglike, args.n_sigma, args.min_loglike)
        for (pf, pl, thr) in peaks:
            match = match_known(pf, known, args.match_tol_hz)
            rows.append({
                "dm":           dm,
                "Nt":           nt,
                "subband_f0":   f0_sub,
                "peak_freq_hz": pf,
                "peak_loglike": pl,
                "threshold":    thr,
                "known_match":  match,
            })

    if n_skipped:
        print(f"  ({n_skipped} files skipped due to read errors or insufficient data)")

    rows.sort(key=lambda r: r["peak_loglike"], reverse=True)

    # ------------------------------------------------------------------
    # Write candidates_raw.csv
    # ------------------------------------------------------------------
    raw_cols = ["dm", "Nt", "subband_f0", "peak_freq_hz",
                "peak_loglike", "threshold", "known_match"]
    raw_csv = os.path.join(args.out_dir, "candidates_raw.csv")
    write_csv(raw_csv, raw_cols, rows)

    n_known_raw = sum(1 for r in rows if r["known_match"])
    print(f"\nWrote {raw_csv}")
    print(f"  total peaks above threshold : {len(rows)}")
    print(f"  matched to known pulsars    : {n_known_raw}")
    print(f"  unmatched                   : {len(rows) - n_known_raw}")

    # ------------------------------------------------------------------
    # Write candidates_dedup.csv
    # ------------------------------------------------------------------
    clusters = dedup_candidates(rows, args.dedup_tol_hz)

    # Column order chosen for readability: identification first, then DM
    # diagnostic columns, then bookkeeping.
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
    print(f"  distinct candidates (tol={args.dedup_tol_hz} Hz) : {len(clusters)}")
    print(f"  distinct known pulsars recovered                  : {len(known_names)}")
    if known_names:
        print(f"    {', '.join(known_names)}")
    print(f"  distinct unmatched candidates                     : {len(clusters) - n_known_d}")

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