#!/usr/bin/env python
"""
collate_fold_snr.py

Post-fold SNR collation for the 47 Tuc Viterbi blind search.

Runs psrstat in batches on all .ar files, parses the output,
cross-matches against candidates_dedup.csv, and writes a CSV sorted by
best S/N descending.

Assumes psrstat is already on PATH (sourced by the calling SLURM script
via setup_viterbi_psr.sh before Python is invoked).

Usage
-----
    python -u collate_fold_snr.py \
        --fold-dir  <stage4_fold/> \
        --dedup-csv <candidates_dedup.csv> \
        --setup     <setup_viterbi_psr.sh> \
        --out-csv   <fold_snr_ranked.csv> \
        [--match-tol-hz 0.05] \
        [--batch-size 50]
"""

import argparse
import csv
import glob
import os
import re
import subprocess
import sys


# ---------------------------------------------------------------------------
# psrstat helpers
# ---------------------------------------------------------------------------

def _parse_psrstat_output(stdout):
    """Parse psrstat -Q output lines into {path: snr} dict."""
    snr_map = {}
    for line in stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        path = parts[0]
        snr = None
        for p in parts[1:]:
            m = re.search(r'snr\s*=\s*([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?)', p)
            if m:
                snr = float(m.group(1))
                break
        if snr is None:
            for p in reversed(parts):
                try:
                    snr = float(p)
                    break
                except ValueError:
                    pass
        snr_map[path] = snr
    return snr_map


def run_psrstat_all(ar_files, batch_size=50):
    """
    Run psrstat -c snr -Q in batches to avoid ARG_MAX shell limit.
    psrstat must already be on PATH.
    Returns dict: ar_path -> snr (float or None).
    """
    if not ar_files:
        return {}

    batches = [ar_files[i:i + batch_size]
               for i in range(0, len(ar_files), batch_size)]

    print(f"Running psrstat on {len(ar_files)} .ar files "
          f"in {len(batches)} batches of {batch_size} ...", flush=True)

    snr_map = {}
    for b_idx, batch in enumerate(batches):
        print(f"  batch {b_idx + 1}/{len(batches)} ({len(batch)} files) ...",
              flush=True)
        file_list = " ".join(f'"{f}"' for f in batch)
        cmd = f'psrstat -c snr -Q {file_list}'
        try:
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
        except subprocess.TimeoutExpired:
            print(f"  WARNING: psrstat timed out on batch {b_idx + 1}.",
                  file=sys.stderr)
            continue

        if result.returncode != 0:
            print(f"  WARNING: psrstat exited {result.returncode} "
                  f"on batch {b_idx + 1}", file=sys.stderr)
            print(f"  stderr: {result.stderr[:300]}", file=sys.stderr)

        snr_map.update(_parse_psrstat_output(result.stdout))

    return snr_map


# ---------------------------------------------------------------------------
# Discover .ar files and group by candidate
# ---------------------------------------------------------------------------

def find_ar_files(fold_dir):
    """
    Returns:
        all_ar:   sorted list of all .ar paths
        cand_map: dict cand_id -> {"poly": [paths], "kepler": [paths]}

    Handles both old layout (cand_NNN/*.ar) and new layout
    (cand_NNN/poly/*.ar, cand_NNN/kepler/*.ar).
    """
    pattern = os.path.join(fold_dir, "**", "*.ar")
    all_ar = sorted(glob.glob(pattern, recursive=True))

    cand_map = {}
    for ar in all_ar:
        parent      = os.path.basename(os.path.dirname(ar))
        grandparent = os.path.basename(
                          os.path.dirname(os.path.dirname(ar)))

        # New layout: cand_NNN/poly/*.ar or cand_NNN/kepler/*.ar
        if grandparent.startswith("cand_") and parent in ("poly", "kepler"):
            cand_id = grandparent
            label   = parent
        # Old layout: cand_NNN/*.ar (label from filename)
        elif parent.startswith("cand_"):
            cand_id = parent
            fname   = os.path.basename(ar)
            if "_poly" in fname:
                label = "poly"
            elif "_kepler" in fname:
                label = "kepler"
            else:
                label = "poly"
        else:
            continue

        if cand_id not in cand_map:
            cand_map[cand_id] = {"poly": [], "kepler": []}
        cand_map[cand_id][label].append(ar)

    return all_ar, cand_map


# ---------------------------------------------------------------------------
# Parse F0 from candfile
# ---------------------------------------------------------------------------

def parse_f0_from_candfile(cand_dir):
    pattern = os.path.join(cand_dir, "*_psrfold.candfile")
    files = glob.glob(pattern)
    if not files:
        return None, None
    rows = []
    with open(files[0]) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 4:
                try:
                    rows.append(float(parts[3]))
                except ValueError:
                    pass
    f0_poly   = rows[0] if len(rows) > 0 else None
    f0_kepler = rows[1] if len(rows) > 1 else None
    return f0_poly, f0_kepler


# ---------------------------------------------------------------------------
# Load and cross-match dedup CSV
# ---------------------------------------------------------------------------

def load_dedup(dedup_csv):
    rows = []
    with open(dedup_csv) as fh:
        lines = [l.rstrip("\n") for l in fh]
    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.split(",")
        if len(parts) < 7:
            continue
        rows.append({
            "peak_freq_hz": float(parts[0]),
            "peak_loglike": float(parts[1]),
            "multiplicity": parts[2],
            "nt_values":    ",".join(parts[3:-4]),
            "Nt_best":      parts[-4],
            "subband_f0":   parts[-3],
            "threshold":    parts[-2],
            "known_match":  parts[-1].strip(),
        })
    return rows


def match_dedup(freq_hz, dedup_rows, tol_hz):
    best, best_df = None, tol_hz
    for row in dedup_rows:
        df = abs(row["peak_freq_hz"] - freq_hz)
        if df < best_df:
            best, best_df = row, df
    return best


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fold-dir",     required=True,
                    help="Root of stage4_fold output (cand_NNN/ subdirs with .ar files).")
    ap.add_argument("--cand-dir",     required=True,
                    help="Root of stage3_viterbi/candidates/ (cand_NNN/ subdirs with candfiles).")
    ap.add_argument("--dedup-csv",    required=True)
    ap.add_argument("--setup",        required=True,
                    help="Path to setup script (must be sourced by caller).")
    ap.add_argument("--out-csv",      required=True)
    ap.add_argument("--match-tol-hz", type=float, default=0.05)
    ap.add_argument("--batch-size",   type=int,   default=50)
    args = ap.parse_args()

    for p, label in [(args.fold_dir,  "fold-dir"),
                     (args.cand_dir,  "cand-dir"),
                     (args.dedup_csv, "dedup-csv"),
                     (args.setup,     "setup")]:
        if not os.path.exists(p):
            sys.exit(f"ERROR: {label} not found: {p}")

    print(f"Loading dedup CSV: {args.dedup_csv}", flush=True)
    dedup_rows = load_dedup(args.dedup_csv)
    print(f"  {len(dedup_rows)} candidates.", flush=True)

    print(f"\nDiscovering .ar files under {args.fold_dir} ...", flush=True)
    all_ar, cand_map = find_ar_files(args.fold_dir)
    print(f"  {len(all_ar)} .ar files in {len(cand_map)} candidate dirs.",
          flush=True)

    if not all_ar:
        sys.exit("ERROR: no .ar files found.")

    snr_map = run_psrstat_all(all_ar, batch_size=args.batch_size)
    print(f"  psrstat returned S/N for {len(snr_map)} files.", flush=True)

    results = []
    for cand_id, ar_dict in sorted(cand_map.items()):
        # candfiles are in cand_dir (stage3_viterbi/candidates/),
        # not in fold_dir (stage4_fold/)
        candfile_dir = os.path.join(args.cand_dir, cand_id)
        f0_poly, f0_kepler = parse_f0_from_candfile(candfile_dir)

        snr_poly, ar_poly_best = None, None
        for ar in ar_dict["poly"]:
            s = snr_map.get(ar)
            if s is not None and (snr_poly is None or s > snr_poly):
                snr_poly, ar_poly_best = s, ar

        snr_kepler, ar_kepler_best = None, None
        for ar in ar_dict["kepler"]:
            s = snr_map.get(ar)
            if s is not None and (snr_kepler is None or s > snr_kepler):
                snr_kepler, ar_kepler_best = s, ar

        if snr_poly is not None and snr_kepler is not None:
            if snr_poly >= snr_kepler:
                snr_best, best_fit, best_ar = snr_poly,   "poly",   ar_poly_best
            else:
                snr_best, best_fit, best_ar = snr_kepler, "kepler", ar_kepler_best
        elif snr_poly is not None:
            snr_best, best_fit, best_ar = snr_poly,   "poly",   ar_poly_best
        elif snr_kepler is not None:
            snr_best, best_fit, best_ar = snr_kepler, "kepler", ar_kepler_best
        else:
            snr_best, best_fit, best_ar = None, "none", None

        match_freq = f0_poly if f0_poly is not None else f0_kepler
        dm = match_dedup(match_freq, dedup_rows, args.match_tol_hz) \
             if match_freq else None

        results.append({
            "cand_id":      cand_id,
            "freq_hz":      float(dm["peak_freq_hz"]) if dm else (match_freq or ""),
            "snr_best":     snr_best      if snr_best  is not None else "",
            "best_fit":     best_fit,
            "snr_poly":     snr_poly      if snr_poly  is not None else "",
            "snr_kepler":   snr_kepler    if snr_kepler is not None else "",
            "peak_loglike": float(dm["peak_loglike"]) if dm else "",
            "multiplicity": dm["multiplicity"]         if dm else "",
            "nt_values":    dm["nt_values"]            if dm else "",
            "known_match":  dm["known_match"]          if dm else "",
            "best_ar":      best_ar if best_ar else "",
        })

    results.sort(key=lambda r: -float(r["snr_best"])
                 if r["snr_best"] != "" else float("inf"))

    cols = ["cand_id", "freq_hz", "snr_best", "best_fit",
            "snr_poly", "snr_kepler", "peak_loglike",
            "multiplicity", "nt_values", "known_match", "best_ar"]

    with open(args.out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in results:
            w.writerow({c: r[c] for c in cols})

    n_known  = sum(1 for r in results if r["known_match"])
    n_no_snr = sum(1 for r in results if r["snr_best"] == "")
    print(f"\nWrote {args.out_csv}", flush=True)
    print(f"  {len(results)} candidates", flush=True)
    print(f"  known pulsar matches: {n_known}", flush=True)
    print(f"  missing S/N:          {n_no_snr}", flush=True)

    print("\nTop 20 by S/N:", flush=True)
    print(f"  {'cand_id':<12} {'freq_hz':>10} {'snr_best':>9} "
          f"{'best_fit':>8} {'loglike':>9} {'mult':>5}  known_match",
          flush=True)
    for r in results[:20]:
        print(f"  {r['cand_id']:<12} {str(r['freq_hz']):>10} "
              f"{str(r['snr_best']):>9} {r['best_fit']:>8} "
              f"{str(r['peak_loglike']):>9} {str(r['multiplicity']):>5}  "
              f"{r['known_match']}", flush=True)


if __name__ == "__main__":
    main()