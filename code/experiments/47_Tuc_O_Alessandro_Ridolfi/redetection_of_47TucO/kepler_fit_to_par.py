#!/usr/bin/env python3
"""
Build psrfold_fil-ready binary .par files from a Viterbi Kepler fit.

Reads the fit summary (kepler_f_spin_hz, kepler_K_hz, kepler_Pb_s,
kepler_T0_s) and the observation start MJD, and writes BT-model pars.

Because the fit-model T0 (where the sine argument is zero) and the BT
periastron/ascending-node T0 differ by a quarter orbit of unknown sign,
this writes THREE pars: T0 as-is, T0 + Pb/4, T0 - Pb/4. Fold all three;
the correct convention is the one giving the sharp, high-S/N fold.
"""
import argparse
import math
import re


def read_fit_summary(path):
    vals = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            parts = line.split()
            key = parts[0]
            if key in ("kepler_f_spin_hz", "kepler_K_hz",
                       "kepler_Pb_s", "kepler_T0_s"):
                vals[key] = float(parts[1])
    missing = {"kepler_f_spin_hz", "kepler_K_hz",
               "kepler_Pb_s", "kepler_T0_s"} - set(vals)
    if missing:
        raise ValueError(f"Missing keys in fit summary: {missing}")
    return vals


def read_mjd_start_from_inf(path):
    with open(path) as fh:
        for line in fh:
            if "Epoch of observation" in line:
                return float(line.split("=")[1].strip())
    raise ValueError("Could not find 'Epoch of observation' in .inf")


PAR_TEMPLATE = """PSRJ            {srcname}
RAJ             {raj}
DECJ            {decj}
F0              {f0:.10f}
DM              {dm}
PEPOCH          {t0:.11f}
BINARY          BT
A1              {a1:.7f}
E               0.0
OM              0.0
T0              {t0:.11f}
PB              {pb_days:.12f}
CLK             UTC(NIST)
EPHEM           DE421
UNITS           TDB
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit-summary", required=True,
                    help="Path to *_fit_summary.dat")
    ap.add_argument("--inf", required=True,
                    help="Path to the PRESTO .inf (for mjd_start)")
    ap.add_argument("--raj", default="00:24:04.65254220")
    ap.add_argument("--decj", default="-72:04:53.7669676")
    ap.add_argument("--dm", default="24.356441")
    ap.add_argument("--srcname", default="J0024-7204O_vitkep")
    ap.add_argument("--out-prefix", required=True,
                    help="Output par prefix; writes <prefix>_asis.par etc.")
    args = ap.parse_args()

    fit = read_fit_summary(args.fit_summary)
    mjd_start = read_mjd_start_from_inf(args.inf)

    f_spin = fit["kepler_f_spin_hz"]
    K = fit["kepler_K_hz"]
    pb_s = fit["kepler_Pb_s"]
    t0_s = fit["kepler_T0_s"]

    a1 = (K / f_spin) * pb_s / (2.0 * math.pi)   # lt-s
    pb_days = pb_s / 86400.0
    t0_mjd_asis = mjd_start + t0_s / 86400.0

    shifts = {
        "asis": 0.0,
        "plus_qtr": pb_s / 4.0,
        "minus_qtr": -pb_s / 4.0,
    }

    print(f"mjd_start   = {mjd_start:.11f}")
    print(f"f_spin (Hz) = {f_spin:.10f}")
    print(f"K (Hz)      = {K:.6e}")
    print(f"Pb (s)      = {pb_s:.6f}  ({pb_days:.12f} d)")
    print(f"a1 sin i    = {a1:.7f} lt-s")
    print(f"T0 as-is    = {t0_mjd_asis:.11f} MJD")
    print()

    for label, shift_s in shifts.items():
        t0 = mjd_start + (t0_s + shift_s) / 86400.0
        par = PAR_TEMPLATE.format(
            srcname="J0024-7204O",
            raj=args.raj, decj=args.decj,
            f0=f_spin, dm=args.dm,
            a1=a1, t0=t0, pb_days=pb_days,
        )
        out = f"{args.out_prefix}_{label}.par"
        with open(out, "w") as fh:
            fh.write(par)
        print(f"wrote {out}   (T0 = {t0:.11f}, shift {shift_s:+.3f} s)")


if __name__ == "__main__":
    main()