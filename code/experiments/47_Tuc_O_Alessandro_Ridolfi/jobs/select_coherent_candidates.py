#!/usr/bin/env python
"""
select_coherent_candidates.py

Builds a shortlist of candidates worth visual inspection using two criteria:

1. BINARY candidates: path_std_hz > binary_threshold (default 0.01 Hz).
   These have a sinusoidal Viterbi path indicating orbital Doppler motion.

2. ISOLATED candidates: peak_loglike in the top N (default top 30) among
   candidates that did NOT pass the binary threshold.
   These have flat paths but high log-likelihood.

The two sets are merged, duplicates removed, and sorted by snr_best descending.
A full _all.csv with path_std for every candidate is also written.

Usage
-----
    python select_coherent_candidates.py \
        --cand-dir        <stage3_viterbi/candidates/> \
        --snr-csv         <stage4_fold/fold_snr_ranked.csv> \
        --out-csv         <shortlist_coherent.csv> \
        [--binary-threshold-hz  0.01] \
        [--isolated-top-n       30]
"""

import argparse
import csv
import glob
import os
import sys
import numpy as np


def read_track(track_path):
    """Return Viterbi path frequencies from _track.dat, or None on failure."""
    try:
        data = np.loadtxt(track_path, comments='#')
        if data.ndim < 2 or data.shape[1] < 2:
            return None
        return data[:, 1]
    except Exception as exc:
        print(f"  WARNING: could not read {track_path}: {exc}", file=sys.stderr)
        return None


def load_snr_csv(snr_csv_path):
    """Load fold_snr_ranked.csv into a dict keyed by cand_id."""
    rows = {}
    if not os.path.isfile(snr_csv_path):
        print(f"WARNING: SNR CSV not found: {snr_csv_path}", file=sys.stderr)
        return rows
    with open(snr_csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows[row["cand_id"]] = row
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cand-dir",             required=True)
    ap.add_argument("--snr-csv",              required=True)
    ap.add_argument("--out-csv",              required=True)
    ap.add_argument("--binary-threshold-hz",  type=float, default=0.01,
                    help="Minimum path std (Hz) to flag as binary (default 0.01).")
    ap.add_argument("--isolated-top-n",       type=int,   default=30,
                    help="Number of top-loglike flat-path candidates to include "
                         "(default 30).")
    args = ap.parse_args()

    if not os.path.isdir(args.cand_dir):
        sys.exit(f"ERROR: cand-dir not found: {args.cand_dir}")

    print(f"Loading SNR CSV: {args.snr_csv}", flush=True)
    snr_map = load_snr_csv(args.snr_csv)
    print(f"  {len(snr_map)} entries.", flush=True)

    cand_dirs = sorted(glob.glob(os.path.join(args.cand_dir, "cand_*")))
    print(f"\nFound {len(cand_dirs)} candidate directories.", flush=True)

    results = []

    for cand_dir in cand_dirs:
        cand_id = os.path.basename(cand_dir)

        track_files = glob.glob(os.path.join(cand_dir, "*_track.dat"))
        if not track_files:
            continue

        path_freqs = read_track(track_files[0])
        if path_freqs is None or len(path_freqs) < 3:
            continue

        path_std   = float(np.std(path_freqs))
        path_mean  = float(np.mean(path_freqs))
        path_range = float(np.max(path_freqs) - np.min(path_freqs))

        snr_row = snr_map.get(cand_id, {})

        # Parse peak_loglike safely
        try:
            peak_loglike = float(snr_row.get("peak_loglike", "nan"))
        except ValueError:
            peak_loglike = float("nan")

        try:
            snr_best = float(snr_row.get("snr_best", "nan"))
        except ValueError:
            snr_best = float("nan")

        results.append({
            "cand_id":       cand_id,
            "path_std_hz":   round(path_std,   6),
            "path_mean_hz":  round(path_mean,  6),
            "path_range_hz": round(path_range, 6),
            "peak_loglike":  peak_loglike,
            "freq_hz":       snr_row.get("freq_hz",      ""),
            "snr_best":      snr_best,
            "best_fit":      snr_row.get("best_fit",     ""),
            "snr_poly":      snr_row.get("snr_poly",     ""),
            "snr_kepler":    snr_row.get("snr_kepler",   ""),
            "multiplicity":  snr_row.get("multiplicity", ""),
            "nt_values":     snr_row.get("nt_values",    ""),
            "known_match":   snr_row.get("known_match",  ""),
            "best_ar":       snr_row.get("best_ar",      ""),
            "track_path":    track_files[0],
        })

    print(f"Processed {len(results)} candidates.", flush=True)

    # --- write full _all.csv ------------------------------------------------
    cols = ["cand_id", "path_std_hz", "path_mean_hz", "path_range_hz",
            "freq_hz", "snr_best", "best_fit", "snr_poly", "snr_kepler",
            "peak_loglike", "multiplicity", "nt_values", "known_match",
            "best_ar", "track_path"]

    all_csv = args.out_csv.replace(".csv", "_all.csv")
    results_sorted_std = sorted(results, key=lambda r: r["path_std_hz"])
    with open(all_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in results_sorted_std:
            w.writerow({c: r[c] for c in cols})
    print(f"\nFull results: {all_csv}", flush=True)

    # --- criterion 1: binary (path_std above threshold) ---------------------
    binary_cands = [r for r in results
                    if r["path_std_hz"] >= args.binary_threshold_hz]
    binary_cands.sort(key=lambda r: r["path_std_hz"], reverse=True)
    print(f"\nBinary candidates (path_std >= {args.binary_threshold_hz} Hz): "
          f"{len(binary_cands)}", flush=True)

    # --- criterion 2: isolated (flat path, top N by peak_loglike) -----------
    flat_cands = [r for r in results
                  if r["path_std_hz"] < args.binary_threshold_hz
                  and np.isfinite(r["peak_loglike"])]
    flat_cands.sort(key=lambda r: r["peak_loglike"], reverse=True)
    isolated_cands = flat_cands[:args.isolated_top_n]
    print(f"Isolated candidates (top {args.isolated_top_n} by loglike "
          f"among flat-path): {len(isolated_cands)}", flush=True)

    # --- merge and sort by snr_best -----------------------------------------
    seen = set()
    shortlist = []
    for r in binary_cands + isolated_cands:
        if r["cand_id"] not in seen:
            seen.add(r["cand_id"])
            shortlist.append(r)

    shortlist.sort(key=lambda r: -r["snr_best"] if np.isfinite(r["snr_best"])
                   else float("inf"))

    print(f"\nTotal shortlist: {len(shortlist)} candidates", flush=True)

    # Add a label column
    binary_ids  = {r["cand_id"] for r in binary_cands}
    isolated_ids = {r["cand_id"] for r in isolated_cands}
    for r in shortlist:
        if r["cand_id"] in binary_ids and r["cand_id"] in isolated_ids:
            r["selection"] = "binary+isolated"
        elif r["cand_id"] in binary_ids:
            r["selection"] = "binary"
        else:
            r["selection"] = "isolated"

    # --- write shortlist CSV ------------------------------------------------
    cols_short = ["cand_id", "selection", "path_std_hz", "path_range_hz",
                  "freq_hz", "snr_best", "best_fit", "snr_poly", "snr_kepler",
                  "peak_loglike", "multiplicity", "nt_values", "known_match",
                  "best_ar", "track_path"]

    with open(args.out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols_short)
        w.writeheader()
        for r in shortlist:
            w.writerow({c: r.get(c, "") for c in cols_short})
    print(f"Shortlist CSV: {args.out_csv}", flush=True)

    # --- summary print ------------------------------------------------------
    print(f"\nTop 30 shortlist (sorted by S/N):", flush=True)
    print(f"  {'cand_id':<12} {'select':<10} {'path_std':>9} "
          f"{'freq_hz':>10} {'snr_best':>9} {'loglike':>9} known_match",
          flush=True)
    for r in shortlist[:30]:
        print(f"  {r['cand_id']:<12} "
              f"{r['selection']:<10} "
              f"{r['path_std_hz']:>9.4f} "
              f"{str(r['freq_hz']):>10} "
              f"{str(r['snr_best']):>9} "
              f"{str(r['peak_loglike']):>9} "
              f"{r['known_match']}", flush=True)

    # --- std distribution ---------------------------------------------------
    stds = [r["path_std_hz"] for r in results]
    print(f"\nPath std distribution:", flush=True)
    for thr in [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]:
        n = sum(1 for s in stds if s >= thr)
        print(f"  >= {thr:.3f} Hz: {n} candidates", flush=True)


if __name__ == "__main__":
    main()