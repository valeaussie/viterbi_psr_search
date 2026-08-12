#!/usr/bin/env python3
"""Merge the per-N_T threshold rows from the array job into one table.

Reads Lth_grid_parts/Lth_Nt*.csv (written by compute_Lth.py --out), reads
L_peak directly from each run's *_loglikes.dat, and writes Lth_grid.csv with
the reportable detection surface det_surf = L_peak - L_th(N_T). Fold S/N is
merged from norm_viterbi_significance.csv if that file is present.
"""
import csv
import glob
import os
import argparse
import numpy as np

PREFIX = "47TucO_redetect_Nt"
SUFFIX = "_loglikes.dat"

def read_fold_snr(path):
    d = {}
    if path and os.path.exists(path):
        for r in csv.DictReader(open(path)):
            try:
                d[int(r["N_T"])] = float(r["fold_snr"])
            except (ValueError, KeyError):
                pass
    return d

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parts", default="Lth_grid_parts/Lth_Nt*.csv")
    ap.add_argument("--dir", default=".", help="dir holding the *_loglikes.dat")
    ap.add_argument("--norm-csv", default="norm_viterbi_significance.csv")
    ap.add_argument("--out", default="Lth_grid.csv")
    args = ap.parse_args()

    fold = read_fold_snr(args.norm_csv)
    rows = []
    for p in sorted(glob.glob(args.parts)):
        r = next(csv.DictReader(open(p)))
        NT = int(r["N_T"])
        llfile = os.path.join(args.dir, f"{PREFIX}{NT}{SUFFIX}")
        Lpeak = float(np.loadtxt(llfile)[:, 1].max())
        Lth = float(r["L_th"])
        rows.append({
            "N_T": NT,
            "T_coh_s": float(r["T_coh_s"]),
            "Nbin": int(r["Nbin"]),
            "L_peak": round(Lpeak, 3),
            "L_th": Lth,
            "sigma_Lth": float(r["sigma_Lth"]),
            "det_surf": round(Lpeak - Lth, 3),
            "fold_snr": fold.get(NT, ""),
        })
    rows.sort(key=lambda x: x["N_T"])

    cols = ["N_T", "T_coh_s", "Nbin", "L_peak", "L_th", "sigma_Lth", "det_surf", "fold_snr"]
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print(" ".join(f"{c:>10}" for c in cols))
    for r in rows:
        print(" ".join(f"{str(r[c]):>10}" for c in cols))
    print(f"\nwrote {args.out}")

if __name__ == "__main__":
    main()