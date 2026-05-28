#!/usr/bin/env python
"""
aggregate_blind.py

Aggregate the blind Viterbi search results. Reads every
*_loglike_curve.dat file under the blind output tree, finds peaks in each
loglike-vs-frequency curve, assembles a master candidate table, and
cross-matches candidate frequencies against the known 47 Tuc pulsar
catalogue (from Stage 0).

A peak is any local maximum in the loglike curve that exceeds a threshold.
For this first pass the threshold is a simple robust cut:

    L_peak > median(L) + n_sigma * 1.4826 * MAD(L)

computed per curve (per subband, per Nt). This is a placeholder until the
formal L_th from Stage 4 is available; it is deliberately conservative and
only used to build the candidate list for inspection.

Usage
-----
    python aggregate_blind.py \
        --blind-dir <stage3_viterbi/blind_v1> \
        --known-yaml <config/known_47tuc_pulsars.yaml> \
        --out-csv <candidates.csv> \
        [--n-sigma 8] [--min-loglike 0] [--match-tol-hz 0.05]
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


def find_peaks(freqs, loglike, n_sigma, min_loglike, min_separation_bins=5):
    """
    Find local maxima in a loglike curve exceeding a robust threshold.

    Returns a list of (peak_freq, peak_loglike, threshold) tuples.
    """
    L = np.asarray(loglike, dtype=float)
    f = np.asarray(freqs, dtype=float)
    if L.size < 3:
        return []

    med = np.median(L)
    mad = np.median(np.abs(L - med))
    sigma = 1.4826 * mad if mad > 0 else np.std(L)
    thr = med + n_sigma * sigma
    thr = max(thr, min_loglike)

    peaks = []
    last_peak_idx = -10**9
    # local maxima above threshold
    for i in range(1, L.size - 1):
        if L[i] >= L[i - 1] and L[i] > L[i + 1] and L[i] > thr:
            if i - last_peak_idx >= min_separation_bins:
                peaks.append((float(f[i]), float(L[i]), float(thr)))
                last_peak_idx = i
            else:
                # keep the higher of two close peaks
                if peaks and L[i] > peaks[-1][1]:
                    peaks[-1] = (float(f[i]), float(L[i]), float(thr))
                    last_peak_idx = i
    return peaks


def parse_run_metadata(path):
    """
    Extract Nt and f0 from a path like .../Nt32/f0_378.000/blind_Nt32_f0_378.000_loglike_curve.dat
    """
    nt = None
    f0 = None
    m = re.search(r"/Nt(\d+)/", path)
    if m:
        nt = int(m.group(1))
    m = re.search(r"/f0_([0-9.]+)/", path)
    if m:
        f0 = float(m.group(1))
    return nt, f0


def load_known_pulsars(yaml_path):
    """Return a list of (name, F0_hz) for known pulsars with an F0."""
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
    """Return the name of a known pulsar within tol_hz of freq, else None."""
    best = None
    best_df = tol_hz
    for name, f0 in known:
        df = abs(freq - f0)
        if df <= best_df:
            best = name
            best_df = df
    return best


def dedup_candidates(rows, tol_hz):
    """
    Collapse peaks within tol_hz of each other into single distinct
    candidates, merging across Nt values and overlapping subbands.

    Greedy: process peaks in descending loglike order; each peak either
    joins an existing cluster (if within tol_hz of its representative
    frequency) or starts a new one. The representative keeps the max-loglike
    peak's frequency and loglike, and records:
      - multiplicity: how many raw peaks merged
      - nt_values: sorted list of distinct Nt at which it appeared
      - known_match: inherited from any merged peak that matched

    Returns a list of dicts sorted by loglike descending.
    """
    ordered = sorted(rows, key=lambda r: r["peak_loglike"], reverse=True)
    clusters = []  # each: dict with rep info + members

    for r in ordered:
        f = r["peak_freq_hz"]
        placed = False
        for c in clusters:
            if abs(f - c["peak_freq_hz"]) <= tol_hz:
                c["multiplicity"] += 1
                c["nt_set"].add(r["Nt"])
                if not c["known_match"] and r["known_match"]:
                    c["known_match"] = r["known_match"]
                placed = True
                break
        if not placed:
            clusters.append({
                "peak_freq_hz": f,
                "peak_loglike": r["peak_loglike"],
                "Nt_best": r["Nt"],
                "subband_f0": r["subband_f0"],
                "threshold": r["threshold"],
                "known_match": r["known_match"],
                "multiplicity": 1,
                "nt_set": {r["Nt"]},
            })

    for c in clusters:
        c["nt_values"] = ",".join(str(x) for x in sorted(
            v for v in c["nt_set"] if v is not None))
        del c["nt_set"]

    clusters.sort(key=lambda c: c["peak_loglike"], reverse=True)
    return clusters


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blind-dir", required=True,
                    help="Root of the blind output tree (contains Nt*/f0_*/).")
    ap.add_argument("--known-yaml", default=None,
                    help="Path to known_47tuc_pulsars.yaml for cross-matching.")
    ap.add_argument("--out-csv", required=True,
                    help="Output CSV path for the candidate table.")
    ap.add_argument("--n-sigma", type=float, default=8.0,
                    help="Robust peak threshold in MAD-sigma above the median.")
    ap.add_argument("--min-loglike", type=float, default=0.0,
                    help="Absolute minimum loglike for a peak.")
    ap.add_argument("--match-tol-hz", type=float, default=0.05,
                    help="Frequency tolerance for known-pulsar cross-match.")
    ap.add_argument("--dedup-tol-hz", type=float, default=0.5,
                    help="Frequency tolerance for collapsing duplicate peaks "
                         "(across Nt and overlapping subbands) into one "
                         "distinct candidate.")
    args = ap.parse_args()

    pattern = os.path.join(args.blind_dir, "**", "*_loglike_curve.dat")
    files = sorted(glob.glob(pattern, recursive=True))
    if not files:
        print(f"No *_loglike_curve.dat files found under {args.blind_dir}")
        return

    known = load_known_pulsars(args.known_yaml)
    print(f"Loaded {len(known)} known pulsars with F0.")
    print(f"Scanning {len(files)} loglike curves...")

    rows = []
    for path in files:
        nt, f0sub = parse_run_metadata(path)
        try:
            arr = np.loadtxt(path)
        except Exception as exc:
            print(f"  skip {path}: {exc}")
            continue
        if arr.ndim != 2 or arr.shape[0] < 3:
            continue
        freqs, loglike = arr[:, 0], arr[:, 1]
        peaks = find_peaks(freqs, loglike, args.n_sigma, args.min_loglike)
        for (pf, pl, thr) in peaks:
            match = match_known(pf, known, args.match_tol_hz)
            rows.append({
                "Nt": nt,
                "subband_f0": f0sub,
                "peak_freq_hz": pf,
                "peak_loglike": pl,
                "threshold": thr,
                "known_match": match if match else "",
            })

    # sort by loglike descending
    rows.sort(key=lambda r: r["peak_loglike"], reverse=True)

    cols = ["Nt", "subband_f0", "peak_freq_hz", "peak_loglike",
            "threshold", "known_match"]
    with open(args.out_csv, "w") as fh:
        fh.write(",".join(cols) + "\n")
        for r in rows:
            fh.write(",".join(str(r[c]) for c in cols) + "\n")

    n_known = sum(1 for r in rows if r["known_match"])
    n_new = len(rows) - n_known
    print(f"\nWrote {args.out_csv}")
    print(f"  total peaks above threshold: {len(rows)}")
    print(f"  matched to known pulsars:    {n_known}")
    print(f"  unmatched (candidates):      {n_new}")

    # ------------------------------------------------------------------
    # De-duplicated candidate list
    # ------------------------------------------------------------------
    clusters = dedup_candidates(rows, args.dedup_tol_hz)
    dedup_csv = os.path.splitext(args.out_csv)[0] + "_dedup.csv"
    dcols = ["peak_freq_hz", "peak_loglike", "multiplicity", "nt_values",
             "Nt_best", "subband_f0", "threshold", "known_match"]
    with open(dedup_csv, "w") as fh:
        fh.write(",".join(dcols) + "\n")
        for c in clusters:
            fh.write(",".join(str(c[k]) for k in dcols) + "\n")

    n_known_d = sum(1 for c in clusters if c["known_match"])
    n_new_d = len(clusters) - n_known_d
    known_names = sorted({c["known_match"] for c in clusters if c["known_match"]})
    print(f"\nWrote {dedup_csv}")
    print(f"  distinct candidates (dedup, tol={args.dedup_tol_hz} Hz): {len(clusters)}")
    print(f"  distinct known pulsars recovered: {len(known_names)}")
    print(f"  distinct unmatched candidates:    {n_new_d}")

    print("\nTop 20 distinct candidates:")
    print(f"  {'freq_hz':>12} {'loglike':>9} {'mult':>5} {'nt_values':>16}  known_match")
    for c in clusters[:20]:
        print(f"  {c['peak_freq_hz']:>12.4f} {c['peak_loglike']:>9.1f} "
              f"{c['multiplicity']:>5} {c['nt_values']:>16}  {c['known_match']}")


if __name__ == "__main__":
    main()