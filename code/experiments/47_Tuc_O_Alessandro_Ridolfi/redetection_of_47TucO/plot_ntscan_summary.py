#!/usr/bin/env python3
"""Two-panel summary of the coherent-timescale scan of 47 Tuc O.

Reads Lth_grid.csv (from collate_Lth.py). Top panel: the detection statistic
L and the threshold L_th versus T_coh (log y), both rising with N_T while L
stays above L_th at every point. Bottom panel: the detection surface
L - L_th and the fold S/N. Chandler bound marked on both. No title; PDF out.
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
OUT_PDF = "ntscan_summary.pdf"
BOUND_S = 452.0
USETEX  = True
C_L     = "#0f6e56"   # detection statistic L
C_LTH   = "#a32d2d"   # threshold L_th
C_DET   = "#185fa5"   # detection surface L - L_th
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
    d = {k: [] for k in ["T", "L", "Lth", "det", "fold"]}
    with open(path) as f:
        for r in csv.DictReader(f):
            d["T"].append(float(r["T_coh_s"]))
            d["L"].append(float(r["L_peak"]))
            d["Lth"].append(float(r["L_th"]))
            d["det"].append(float(r["det_surf"]))
            d["fold"].append(float(r["fold_snr"]))
    idx = np.argsort(d["T"])
    return {k: np.asarray(v)[idx] for k, v in d.items()}

def L(s, tex, alt):
    return s if tex else alt

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-csv", default=IN_CSV)
    ap.add_argument("--out", default=OUT_PDF)
    ap.add_argument("--bound", type=float, default=BOUND_S)
    ap.add_argument("--no-usetex", action="store_true")
    args = ap.parse_args()
    tex = USETEX and not args.no_usetex
    set_style(tex)
    d = load(args.in_csv)

    fig, (axT, axB) = plt.subplots(2, 1, figsize=(8, 8), sharex=True,
                                   gridspec_kw={"height_ratios": [1, 1], "hspace": 0.08})

    # --- top: L and L_th ---
    axT.set_yscale("log")
    axT.fill_between(d["T"], d["Lth"], d["L"], color=C_L, alpha=0.08)
    axT.plot(d["T"], d["L"], "-o", color=C_L, lw=1.7, ms=6,
             label=L(r"$L$ (recovered)", tex, "L (recovered)"))
    axT.plot(d["T"], d["Lth"], "--s", color=C_LTH, lw=1.5, ms=5,
             label=L(r"$L_{\mathrm{th}}$ (threshold, $\alpha'=0.1$)", tex, "L_th (threshold)"))
    axT.axvline(args.bound, color="0.35", ls=":", lw=1.3)
    axT.set_ylabel(L(r"Log-likelihood", tex, "Log-likelihood"))
    axT.legend(loc="upper left", fontsize=11, framealpha=0.9)

    # --- bottom: det surface (left) and fold S/N (right) ---
    axB.set_xscale("log")
    axB.plot(d["T"], d["det"], "-o", color=C_DET, lw=1.7, ms=6,
             label=L(r"$L-L_{\mathrm{th}}$", tex, "L - L_th"))
    axB.axvline(args.bound, color="0.35", ls=":", lw=1.3)
    axB.text(args.bound * 1.04, axB.get_ylim()[0] + 30,
             L(r"$T_{\mathrm{coh,max}}$", tex, "T_coh,max"),
             color="0.35", rotation=90, va="bottom", fontsize=11)
    axB.set_xlabel(L(r"$T_{\mathrm{coh}}$ (s)", tex, "T_coh (s)"))
    axB.set_ylabel(L(r"Detection surface $L-L_{\mathrm{th}}$", tex, "L - L_th"), color=C_DET)
    axB.tick_params(axis="y", labelcolor=C_DET)
    axB.set_ylim(0, 1.15 * np.nanmax(d["det"]))

    axB2 = axB.twinx()
    axB2.plot(d["T"], d["fold"], "-^", color=C_FOLD, lw=1.4, ms=6, label="Fold S/N")
    axB2.set_ylabel("Fold S/N", color=C_FOLD)
    axB2.tick_params(axis="y", labelcolor=C_FOLD)
    axB2.set_ylim(0, 1.15 * np.nanmax(d["fold"]))

    h1, l1 = axB.get_legend_handles_labels()
    h2, l2 = axB2.get_legend_handles_labels()
    axB.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=11, framealpha=0.9)

    axB.invert_xaxis()   # long segments (bound side) on the left
    fig.savefig(args.out, bbox_inches="tight")
    print(f"wrote {args.out}")

if __name__ == "__main__":
    main()