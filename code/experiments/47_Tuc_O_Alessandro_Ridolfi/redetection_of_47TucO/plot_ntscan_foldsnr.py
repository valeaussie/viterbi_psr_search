#!/usr/bin/env python3
"""Fold S/N versus T_coh for the coherent-timescale scan of 47 Tuc O.

Reads Lth_grid.csv (from collate_Lth.py) and plots only the fold S/N.
Chandler bound marked with a dotted line (omit with --no-bound).
No title; PDF out.
"""
import csv
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rc

# ----------------------------------------------------------------- config
IN_CSV  = "Lth_grid.csv"
OUT_PDF = "ntscan_foldsnr.pdf"
BOUND_S = 452.0
USETEX  = True
C_FOLD  = "#ba7517"   # fold S/N
# ------------------------------------------------------------------------

def set_style(usetex):
    if usetex:
        matplotlib.rcParams["text.latex.preamble"] += (
            r"\usepackage[dvips]{graphicx}\usepackage{amsmath}\usepackage{amssymb}")
    rc("text", usetex=usetex)
    rc("font", size=14.0)
    rc("font", **{"family": "serif"})
    plt.rcParams["figure.facecolor"] = "white"

def load(path):
    T, fold = [], []
    with open(path) as f:
        for r in csv.DictReader(f):
            T.append(float(r["T_coh_s"]))
            s = r["fold_snr"].strip()
            fold.append(float(s) if s else np.nan)
    T, fold = np.asarray(T), np.asarray(fold)
    idx = np.argsort(T)
    return T[idx], fold[idx]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-csv", default=IN_CSV)
    ap.add_argument("--out", default=OUT_PDF)
    ap.add_argument("--bound", type=float, default=BOUND_S)
    ap.add_argument("--no-bound", action="store_true")
    ap.add_argument("--no-usetex", action="store_true")
    args = ap.parse_args()
    tex = USETEX and not args.no_usetex
    set_style(tex)
    T, fold = load(args.in_csv)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.set_xscale("log")
    ax.plot(T, fold, "-^", color=C_FOLD, lw=1.4, ms=6)
    ax.set_ylim(0, 1.15 * np.nanmax(fold))
    if not args.no_bound:
        ax.axvline(args.bound, color="0.35", ls=":", lw=1.3)
        ax.text(args.bound * 1.04, 0.3,
                r"$T_{\mathrm{coh,max}}$" if tex else "T_coh,max",
                color="0.35", rotation=90, va="bottom", fontsize=11)
    ax.set_xlabel(r"$T_{\mathrm{coh}}$ (s)" if tex else "T_coh (s)")
    ax.set_ylabel("Fold S/N")
    ax.invert_xaxis()   # long segments (bound side) on the left
    fig.savefig(args.out, bbox_inches="tight")
    print(f"wrote {args.out}")

if __name__ == "__main__":
    main()