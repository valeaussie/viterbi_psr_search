#!/usr/bin/env python3
"""
make_polyco_fold.py

Folds the 47 Tuc O redetection with ALL Taylor coefficients from the
Viterbi fit (f0 to fN, here N = 5), by writing them into a TEMPO-format
polyco file that prepfold reads with -psr and -polycos.

Why a polyco file: prepfold's command line accepts only -f, -fd, -fdd,
and moving the expansion from the midpoint to the start of the data
loses thousands of cycles of phase. A polyco block is a phase
polynomial about a midpoint time, so the fit can be passed in exactly,
about its own reference time, with no truncation and no change of frame.

Mapping (exact). The fit gives f(t) = sum_k f_k (t - tref)^k / k!, with t
in seconds from the start of the data (topocentric). The pulse phase is

  phi(t) = sum_k f_k (t - tref)^(k+1) / (k+1)!

A polyco block evaluates, with dt = (t - tref) / 60 in minutes,

  phase(dt) = RPHASE + 60 * dt * F0 + sum_i COEFF_i * dt^i

so we set TMID = epoch + tref / 86400, RPHASE = 0, F0 = f_0,
COEFF_0 = COEFF_1 = 0, and COEFF_(k+1) = f_k * 60^(k+1) / (k+1)! for k >= 1.
The number of coefficients is padded with zeros to a multiple of 3,
because PRESTO reads them three per line.

prepfold evaluates polycos at the raw (topocentric UTC) data times, so
no barycentric conversion is needed; the fit was made on topocentric data.

Checks written to the log:
  - the .inf is not barycentred and its epoch matches the filterbank tstart
  - polyco line 1 is at most 88 characters and line 2 exactly 79 (PRESTO
    reads them into fixed-size buffers, as for TEMPO output)
  - round trip: the polyco file is read back with a Python copy of
    PRESTO's getpoly/phcalc logic, and its phase and frequency are compared
    with the Taylor series evaluated directly on 1001 times across the data

Outputs, in --outdir:
  <name>_<fold>.polycos             polyco file for each fold
  <name>_polyco_fold_params.json    every value, with its source
  <name>_polyco_fold_params.log     the same, human-readable
  <commands-file>                   one prepfold command per fold
                                    (--commands-file, default prepfold_polyco_commands.txt)

Usage (from inside folding_with_PRESTO/):
  python3 make_polyco_fold.py
  python3 make_polyco_fold.py --npfact 4 --ndmfact 10 --nsub 64
  python3 make_polyco_fold.py --fit ../47TucO_redetect_Nt16_fit_summary.dat \
      --name 47TucO_Nt16 --outdir ntscan --folds poly --commands-file commands_Nt16.txt
"""

import argparse
import json
import math
import os
import shlex
from datetime import datetime, timezone

from make_prepfold_inputs import (read_fit_summary, read_inf, inf_get,
                                  read_sigproc_header, taylor_coeffs)

HERE = os.path.dirname(os.path.abspath(__file__))
REDETECT = os.path.dirname(HERE)
EXPERIMENT = os.path.dirname(REDETECT)

DEFAULT_FIT = os.path.join(REDETECT, "47TucO_redetect_Nt32_fit_summary.dat")
DEFAULT_INF = os.path.join(REDETECT, "47TucO_redetect_cfbf00000_Plan1_1_DM24.36.inf")
DEFAULT_FIL = os.path.join(EXPERIMENT, "data", "47Tuc_22UL_1of2_L.fil")

FOLDS = {"poly": "poly_", "kepler_taylor": "kepler_"}


def polyco_coeffs(f):
    """COEFF_i for a Taylor series f[k] (k-th derivative of frequency, Hz/s^k)."""
    c = [0.0, 0.0]
    for k in range(1, len(f)):
        c.append(f[k] * 60.0 ** (k + 1) / math.factorial(k + 1))
    while len(c) % 3:
        c.append(0.0)
    if len(c) > 15:
        raise ValueError("PRESTO accepts at most 15 polyco coefficients")
    return c


def write_polyco(path, psrname, imjd, frac_mid, dm, f0, nblk, coeffs, obsfreq):
    """Write one TEMPO-format polyco block. Returns (line 1 length, line 2 length).
    PRESTO reads line 1 into a 90-byte buffer and line 2 into an 80-byte
    buffer, so line 1 must be at most 88 characters and line 2 exactly 79."""
    secs = frac_mid * 86400.0
    hh = int(secs // 3600)
    mm = int((secs - 3600 * hh) // 60)
    ss = secs - 3600 * hh - 60 * mm
    utc = hh * 10000 + mm * 100 + ss
    date = datetime.now(timezone.utc).strftime("%d-%b-%y").upper()
    tmid = f"{imjd:d}" + f"{frac_mid:.15f}"[1:]          # e.g. 59391.211907...
    line1 = f"{psrname} {date} {utc:.2f} {tmid} {dm:.6f} {0.0:.3f} {0.0:.3f}"
    line2 = f"{0.0:20.6f}{f0:18.12f}{0:5d}{nblk:5d}{len(coeffs):5d}{obsfreq:10.3f}" + " " * 16
    with open(path, "w") as fh:
        fh.write(line1 + "\n")
        fh.write(line2 + "\n")
        for i in range(0, len(coeffs), 3):
            fh.write("".join(f"{c:25.17e}" for c in coeffs[i:i + 3]) + "\n")
    return len(line1), len(line2)


def read_polyco_like_presto(path):
    """Python copy of the parts of PRESTO's getpoly() that matter here."""
    with open(path) as fh:
        l1 = fh.readline()
        l2 = fh.readline()
        tok = l1.split()
        tmid_tok = tok[3]
        imjd = int(tmid_tok.split(".")[0])
        fmid = float("0." + tmid_tok.split(".")[1])
        t2 = l2[:79].split()
        rphase, f0 = float(t2[0]), float(t2[1])
        nblk, ncoeff = int(t2[3]), int(t2[4])
        coeffs = []
        for _ in range(ncoeff // 3):
            coeffs += [float(x.replace("D", "e")) for x in fh.readline().split()]
    return imjd, fmid, rphase - math.floor(rphase), f0, nblk, coeffs


def phcalc_like_presto(imjd, fmid, rphase, f0, nblk, coeffs, mjd0, mjd1):
    """Python copy of PRESTO's phcalc(): unwrapped phase (cycles) and frequency."""
    dtmin = ((mjd0 - imjd) + (mjd1 - fmid)) * 1440.0
    if abs(dtmin) >= nblk / 2.0 + 1e-7:
        raise ValueError(f"time {dtmin} min outside polyco span {nblk} min")
    n = len(coeffs)
    freq, phase = 0.0, coeffs[n - 1]
    for i in range(n - 1, 0, -1):
        freq = dtmin * freq + i * coeffs[i]
        phase = dtmin * phase + coeffs[i - 1]
    return phase + rphase + dtmin * 60.0 * f0, f0 + freq / 60.0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fit", default=DEFAULT_FIT)
    ap.add_argument("--inf", default=DEFAULT_INF)
    ap.add_argument("--fil", default=DEFAULT_FIL)
    ap.add_argument("--outdir", default=HERE)
    ap.add_argument("--name", default="47TucO_Nt32")
    ap.add_argument("--psrname", default="J0024-7204O")
    ap.add_argument("--npfact", type=int, default=2)
    ap.add_argument("--ndmfact", type=int, default=3)
    ap.add_argument("--nsub", type=int, default=None)
    ap.add_argument("--ncpus", type=int, default=1)
    ap.add_argument("--no-searchfdd", action="store_true")
    ap.add_argument("--folds", default="poly,kepler_taylor",
                    help="comma-separated subset of: " + ", ".join(FOLDS))
    ap.add_argument("--commands-file", default="prepfold_polyco_commands.txt")
    args = ap.parse_args()
    args.outdir = os.path.abspath(args.outdir)
    folds = [x.strip() for x in args.folds.split(",")]
    for x in folds:
        if x not in FOLDS:
            raise SystemExit(f"ABORT: unknown fold '{x}'; choose from {list(FOLDS)}")

    fit_path = os.path.abspath(args.fit)
    inf_path = os.path.abspath(args.inf)
    fil_path = os.path.abspath(args.fil)
    fit = read_fit_summary(fit_path)
    inf = read_inf(inf_path)
    hdr = read_sigproc_header(fil_path)

    record = {"created_utc": datetime.now(timezone.utc).isoformat(),
              "inputs": {"fit_summary": fit_path, "inf": inf_path, "filterbank": fil_path},
              "values": {}, "checks": {}, "folds": {}}
    V = record["values"]

    def put(key, value, source):
        V[key] = {"value": value, "source": source}
        return value

    # Observation quantities
    epoch_str = put("inf_epoch_mjd_str", inf_get(inf, "Epoch of observation"), inf_path)
    imjd = int(epoch_str.split(".")[0])
    fmjd = float("0." + epoch_str.split(".")[1])
    bary = int(put("inf_barycentred", inf_get(inf, "Barycentered?"), inf_path))
    dm = float(put("inf_dm", inf_get(inf, "Dispersion measure"), inf_path))
    nbins = int(put("inf_nbins", inf_get(inf, "Number of bins"), inf_path))
    dt = float(put("inf_tsamp_s", inf_get(inf, "Width of each time series bin"), inf_path))
    t_obs = put("t_obs_s", nbins * dt, "inf_nbins * inf_tsamp_s")
    flo = float(put("inf_low_chan_mhz", inf_get(inf, "Central freq of low channel"), inf_path))
    nchan = int(put("inf_nchan", inf_get(inf, "Number of channels"), inf_path))
    chw = float(put("inf_chan_bw_mhz", inf_get(inf, "Channel bandwidth"), inf_path))
    obsfreq = put("obsfreq_centre_mhz", flo + 0.5 * (nchan - 1) * chw,
                  "inf_low_chan_mhz + 0.5 * (inf_nchan - 1) * inf_chan_bw_mhz")

    if bary != 0:
        raise SystemExit("ABORT: .inf says the data are barycentred; the fit must be topocentric.")
    tstart = put("fil_tstart_mjd", hdr["tstart"], fil_path)
    diff_s = put("inf_epoch_minus_fil_tstart_s",
                 ((imjd - math.floor(tstart)) + (fmjd - (tstart - math.floor(tstart)))) * 86400.0,
                 "inf epoch minus filterbank tstart")
    record["checks"]["inf_epoch_matches_fil_tstart_within_one_sample"] = abs(diff_s) < dt
    if abs(diff_s) >= dt:
        raise SystemExit(f"ABORT: .inf epoch and filterbank tstart differ by {diff_s} s.")

    os.makedirs(args.outdir, exist_ok=True)
    grid = [t_obs * i / 1000.0 for i in range(1001)]

    for fold in folds:
        prefix = FOLDS[fold]
        f = put(f"{fold}_coeffs_at_tref", taylor_coeffs(fit, prefix), fit_path)
        tref = put(f"{fold}_tref_s", fit[f"{prefix}tref_s"][0], fit_path)
        put(f"{fold}_n_taylor_coeffs", len(f), f"count of {prefix}f0, {prefix}f1, ... in fit summary")

        frac_mid = fmjd + tref / 86400.0
        imid = imjd + int(math.floor(frac_mid))
        frac_mid -= math.floor(frac_mid)
        put(f"{fold}_tmid_mjd", f"{imid:d}" + f"{frac_mid:.15f}"[1:],
            f"inf epoch + {fold}_tref_s / 86400")
        half_span_min = max(tref, t_obs - tref) / 60.0
        nblk = put(f"{fold}_nspan_min", int(math.ceil(2.0 * half_span_min)) + 2,
                   f"ceil(2 * max({fold}_tref_s, t_obs_s - {fold}_tref_s) / 60) + 2")
        coeffs = put(f"{fold}_polyco_coeffs", polyco_coeffs(f),
                     "COEFF_0 = COEFF_1 = 0; COEFF_(k+1) = f_k 60^(k+1)/(k+1)!; zero-padded to a multiple of 3")

        pc_path = os.path.join(args.outdir, f"{args.name}_{fold}.polycos")
        l1len, l2len = write_polyco(pc_path, args.psrname, imid, frac_mid, dm, f[0], nblk, coeffs, obsfreq)
        put(f"{fold}_polyco_file", pc_path, "written by this script")
        record["checks"][f"{fold}_polyco_line1_at_most_88_chars"] = (l1len <= 88)
        record["checks"][f"{fold}_polyco_line2_is_79_chars"] = (l2len == 79)

        # Round trip: read back as PRESTO would and compare with the Taylor series
        pc = read_polyco_like_presto(pc_path)
        ph0, _ = phcalc_like_presto(*pc, imjd, fmjd)
        dphi, df = 0.0, 0.0
        for t in grid:
            ph, fr = phcalc_like_presto(*pc, imjd, fmjd + t / 86400.0)
            tau = t - tref
            phi_direct = sum(f[k] * tau ** (k + 1) / math.factorial(k + 1) for k in range(len(f)))
            phi_direct0 = sum(f[k] * (-tref) ** (k + 1) / math.factorial(k + 1) for k in range(len(f)))
            f_direct = sum(f[k] * tau ** k / math.factorial(k) for k in range(len(f)))
            dphi = max(dphi, abs((ph - ph0) - (phi_direct - phi_direct0)))
            df = max(df, abs(fr - f_direct))
        put(f"{fold}_roundtrip_max_phase_diff_cycles", dphi,
            "polyco read back (Python copy of PRESTO getpoly/phcalc) vs Taylor series, 1001 times over [0, t_obs_s]")
        put(f"{fold}_roundtrip_max_freq_diff_hz", df, "same comparison, spin frequency")
        _, f_start = phcalc_like_presto(*pc, imjd, fmjd)
        put(f"{fold}_polyco_freq_at_start_hz", f_start, "polyco evaluated at t = 0")
        record["checks"][f"{fold}_roundtrip_phase_below_1e-4_cycles"] = dphi < 1e-4

        cmd = ["prepfold", "-filterbank", "-topo", "-noxwin", "-ncpus", str(args.ncpus),
               "-dm", repr(dm), "-npfact", str(args.npfact), "-ndmfact", str(args.ndmfact)]
        if args.nsub is not None:
            cmd += ["-nsub", str(args.nsub)]
        cmd += ["-psr", args.psrname, "-polycos", pc_path]
        if not args.no_searchfdd:
            cmd.append("-searchfdd")
        cmd += ["-o", f"{args.name}_{fold}_polyco", fil_path]
        record["folds"][fold] = shlex.join(cmd)

    stem = os.path.join(args.outdir, f"{args.name}_polyco_fold_params")
    with open(stem + ".json", "w") as fh:
        json.dump(record, fh, indent=2)
    with open(stem + ".log", "w") as fh:
        fh.write(f"# make_polyco_fold.py  {record['created_utc']}\n")
        for k, v in record["inputs"].items():
            fh.write(f"# input {k}: {v}\n")
        for k, v in V.items():
            fh.write(f"{k} = {v['value']!r}    [from: {v['source']}]\n")
        for k, v in record["checks"].items():
            fh.write(f"CHECK {k}: {v}\n")
        for k, v in record["folds"].items():
            fh.write(f"FOLD {k}: {v}\n")
    with open(os.path.join(args.outdir, args.commands_file), "w") as fh:
        for v in record["folds"].values():
            fh.write(v + "\n")
    with open(stem + ".log") as fh:
        print(fh.read())


if __name__ == "__main__":
    main()