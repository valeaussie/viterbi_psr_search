
#!/usr/bin/env python
"""
check_orbit_model.py

Pre-fold consistency check for a binary fold.

Compares three curves over the observation:
  A. Viterbi fit convention:  f = f_spin - K*sin(2*pi*(t - T0_fit)/Pb)
  B. prepfold convention:     f = f_spin - K*cos(2*pi*(t - To_s)/Pb)
     with To_s = (To_MJD - mjd_start)*86400
  C. the actual Viterbi path (top row of the *_paths.dat file)

Prints max |A - B| (should be ~0 if the To mapping is correct) and the
RMS of (C - B) in bin widths (should be of order 1 for a good orbit).

Usage:
    python check_orbit_model.py \
        --paths-file <blind_..._paths.dat> --tsft <Tsft seconds> \
        --f-spin <Hz> --k <Hz> --pb <s> --t0-s <s> \
        --to-mjd <MJD> --inf <matching .inf> [--df-bin <Hz>] \
        [--plot check.png]
"""

import argparse
import math
import sys

import numpy as np

SEC_PER_DAY = 86400.0


def read_mjd_start_from_inf(inf_path):
    with open(inf_path, errors="replace") as fh:
        for line in fh:
            if "epoch of observation" in line.lower():
                return float(line.split("=")[-1].strip())
    raise ValueError(f"No 'Epoch of observation' line in {inf_path}")


def read_top_path(paths_file):
    """
    *_paths.dat rows are: loglike f_1 f_2 ... f_Nsft
    Rows are written in ascending loglike order, so the best path is
    the LAST row. Returns (loglike, np.array of path frequencies).
    """
    best = None
    with open(paths_file) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 2:
                continue
            best = parts   # keep last non-empty row
    if best is None:
        raise ValueError(f"No path rows in {paths_file}")
    return float(best[0]), np.array([float(x) for x in best[1:]])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paths-file", required=True)
    ap.add_argument("--tsft",  type=float, required=True)
    ap.add_argument("--f-spin", type=float, required=True)
    ap.add_argument("--k",      type=float, required=True)
    ap.add_argument("--pb",     type=float, required=True)
    ap.add_argument("--t0-s",   type=float, required=True,
                    help="Fit T0 (kepler_model convention), s from start")
    ap.add_argument("--to-mjd", type=float, required=True,
                    help="prepfold To in MJD (from setup script or grid)")
    ap.add_argument("--inf", required=True,
                    help=".inf matching the .dat the path came from")
    ap.add_argument("--df-bin", type=float, default=None,
                    help="Frequency bin width (Hz); if omitted, estimated "
                         "as 1/Tsft")
    ap.add_argument("--plot", default=None,
                    help="Optional output PNG path")
    args = ap.parse_args()

    mjd_start = read_mjd_start_from_inf(args.inf)
    loglike, path = read_top_path(args.paths_file)
    n = len(path)
    t = args.tsft / 2.0 + np.arange(n) * args.tsft

    w = 2.0 * math.pi / args.pb
    model_fit = args.f_spin - args.k * np.sin(w * (t - args.t0_s))

    to_s = (args.to_mjd - mjd_start) * SEC_PER_DAY
    model_pf = args.f_spin - args.k * np.cos(w * (t - to_s))

    df_bin = args.df_bin if args.df_bin else 1.0 / args.tsft

    dab = float(np.max(np.abs(model_fit - model_pf)))
    rms_pf = float(np.sqrt(np.mean((path - model_pf) ** 2)))
    rms_fit = float(np.sqrt(np.mean((path - model_fit) ** 2)))

    print(f"Path: {n} segments, loglike = {loglike:.2f}, "
          f"Tsft = {args.tsft:.3f} s, df_bin = {df_bin:.6e} Hz")
    print(f"mjd_start (from .inf)          = {mjd_start:.10f}")
    print(f"To as seconds from data start  = {to_s:.3f} s")
    print()
    print(f"max |fit-model - prepfold-model| = {dab:.3e} Hz "
          f"({dab / df_bin:.3f} bins)")
    if dab > 0.1 * df_bin:
        print("  WARNING: conventions disagree at >0.1 bin. The To value")
        print("  does not correspond to the fit T0 via To = T0 + Pb/4.")
        print("  Do not fold until this is resolved.")
    print(f"RMS(path - prepfold-model) = {rms_pf:.3e} Hz "
          f"({rms_pf / df_bin:.2f} bins)")
    print(f"RMS(path - fit-model)      = {rms_fit:.3e} Hz "
          f"({rms_fit / df_bin:.2f} bins)")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(9, 4))
        ax.plot(t, path, "r-", lw=1.5, label="Viterbi path")
        ax.plot(t, model_fit, "b--", lw=1.2, label="Kepler fit convention")
        ax.plot(t, model_pf, "g:", lw=1.8, label="prepfold convention")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Frequency (Hz)")
        ax.legend()
        fig.tight_layout()
        fig.savefig(args.plot, dpi=150)
        print(f"Wrote {args.plot}")


if __name__ == "__main__":
    sys.exit(main())