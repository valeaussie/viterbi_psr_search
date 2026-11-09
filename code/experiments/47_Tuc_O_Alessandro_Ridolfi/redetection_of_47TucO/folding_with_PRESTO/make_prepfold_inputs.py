#!/usr/bin/env python3
"""
make_prepfold_inputs.py

Builds prepfold commands for the 47 Tuc O redetection folds from the
Viterbi fit summary, the topocentric (not barycentred) .inf file and
the raw filterbank header. Every value is computed here at runtime and
written, with its source file, to:

  <name>_fold_params.json   machine-readable record
  <name>_fold_params.log    human-readable record
  prepfold_commands.txt     one prepfold command per line

Folds (all on the raw filterbank, with -topo, because the fit was made
on topocentric data):

  poly           poly_f0..poly_fN, given at tref, re-expanded exactly to
                 t = 0 (prepfold folds about the start of the data),
                 then truncated to -f, -fd, -fdd.
  kepler_taylor  kepler_f0..kepler_fN, same treatment.
  kepler_orbit   circular orbit via -bin, with
                   x   = K * Pb / (2 * pi * f_spin)   projected semi-major axis (lt-s)
                   T_a = T0 + Pb / 4                  ascending node (s from start)
                   To  = epoch + T_a / 86400          passed to -To (MJD)
                 and e = 0, w = 0 (periastron placed at the ascending node).

Time convention: the fit's t is seconds from the start of the .dat, whose
start is the .inf epoch. The script checks that the .inf epoch matches the
filterbank tstart, and that the .inf is not barycentred.

Usage (from inside folding_with_PRESTO/):
  python3 make_prepfold_inputs.py
  python3 make_prepfold_inputs.py --npfact 4 --ndmfact 10 --nsub 64
"""

import argparse
import json
import math
import os
import shlex
import struct
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
REDETECT = os.path.dirname(HERE)
EXPERIMENT = os.path.dirname(REDETECT)

DEFAULT_FIT = os.path.join(REDETECT, "47TucO_redetect_Nt32_fit_summary.dat")
DEFAULT_INF = os.path.join(REDETECT, "47TucO_redetect_cfbf00000_Plan1_1_DM24.36.inf")
DEFAULT_FIL = os.path.join(EXPERIMENT, "data", "47Tuc_22UL_1of2_L.fil")

SIGPROC_TYPES = {
    "telescope_id": "i", "machine_id": "i", "data_type": "i",
    "barycentric": "i", "pulsarcentric": "i", "nbits": "i",
    "nsamples": "i", "nchans": "i", "nifs": "i", "nbeams": "i", "ibeam": "i",
    "tstart": "d", "tsamp": "d", "fch1": "d", "foff": "d", "refdm": "d",
    "az_start": "d", "za_start": "d", "src_raj": "d", "src_dej": "d",
    "period": "d",
    "source_name": "s", "rawdatafile": "s",
}


def read_sigproc_header(path):
    """Read the SIGPROC filterbank header into a dict."""
    hdr = {}
    with open(path, "rb") as fh:
        def read_string():
            (n,) = struct.unpack("<i", fh.read(4))
            return fh.read(n).decode("ascii", errors="replace")
        if read_string() != "HEADER_START":
            raise ValueError(f"{path}: no HEADER_START, not a SIGPROC filterbank")
        while True:
            key = read_string()
            if key == "HEADER_END":
                break
            kind = SIGPROC_TYPES.get(key)
            if kind is None:
                raise ValueError(f"{path}: unknown header key '{key}'")
            if kind == "s":
                hdr[key] = read_string()
            else:
                fmt = "<" + kind
                (hdr[key],) = struct.unpack(fmt, fh.read(struct.calcsize(fmt)))
    return hdr


def read_fit_summary(path):
    """Return {key: (value, error or None)} for poly_* and kepler_* lines."""
    vals = {}
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 2 or parts[0].startswith("#"):
                continue
            key = parts[0]
            if not (key.startswith("poly_") or key.startswith("kepler_")):
                continue
            try:
                value = float(parts[1])
            except ValueError:
                continue
            err = float(parts[2]) if len(parts) > 2 else None
            vals[key] = (value, err)
    return vals


def read_inf(path):
    """Return {label: value string} for 'label = value' lines of a PRESTO .inf."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if "=" in line:
                label, value = line.rsplit("=", 1)
                out[label.strip()] = value.strip()
    return out


def inf_get(inf, prefix):
    matches = [v for k, v in inf.items() if k.startswith(prefix)]
    if len(matches) != 1:
        raise KeyError(f".inf: expected one line starting '{prefix}', found {len(matches)}")
    return matches[0]


def taylor_coeffs(fit, prefix):
    """Collect <prefix>f0, <prefix>f1, ... in order (k-th derivatives at tref)."""
    coeffs = []
    k = 0
    while f"{prefix}f{k}" in fit:
        coeffs.append(fit[f"{prefix}f{k}"][0])
        k += 1
    if len(coeffs) < 3:
        raise ValueError(f"fewer than 3 coefficients found for prefix '{prefix}'")
    return coeffs


def reexpand(coeffs, delta):
    """Exact re-expansion of a Taylor series about tref to tref + delta.
    coeffs[k] is the k-th derivative at tref; returns k-th derivatives at
    tref + delta."""
    n = len(coeffs)
    return [sum(coeffs[k] * delta ** (k - j) / math.factorial(k - j)
                for k in range(j, n)) for j in range(n)]


def dropped_phase(coeffs_at_start, t):
    """Pulse phase (cycles) at time t from the terms prepfold cannot take
    (derivative order 3 and above), with coefficients referenced to t = 0."""
    return sum(coeffs_at_start[k] * t ** (k + 1) / math.factorial(k + 1)
               for k in range(3, len(coeffs_at_start)))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fit", default=DEFAULT_FIT)
    ap.add_argument("--inf", default=DEFAULT_INF)
    ap.add_argument("--fil", default=DEFAULT_FIL)
    ap.add_argument("--outdir", default=HERE)
    ap.add_argument("--name", default="47TucO_Nt32")
    ap.add_argument("--npfact", type=int, default=2)
    ap.add_argument("--ndmfact", type=int, default=3)
    ap.add_argument("--nsub", type=int, default=None)
    ap.add_argument("--ncpus", type=int, default=1)
    ap.add_argument("--no-searchfdd", action="store_true",
                    help="do not search fdd in the Taylor folds")
    args = ap.parse_args()

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

    # Observation quantities from the .inf
    epoch_str = put("inf_epoch_mjd_str", inf_get(inf, "Epoch of observation"), inf_path)
    imjd_str, frac_str = epoch_str.split(".")
    imjd, fmjd = int(imjd_str), float("0." + frac_str)
    bary = int(put("inf_barycentred", inf_get(inf, "Barycentered?"), inf_path))
    dm = float(put("inf_dm", inf_get(inf, "Dispersion measure"), inf_path))
    nbins = int(put("inf_nbins", inf_get(inf, "Number of bins"), inf_path))
    dt = float(put("inf_tsamp_s", inf_get(inf, "Width of each time series bin"), inf_path))
    t_obs = put("t_obs_s", nbins * dt, "inf_nbins * inf_tsamp_s")

    # Checks
    if bary != 0:
        raise SystemExit("ABORT: .inf says the data are barycentred, but these folds use -topo.")
    tstart = put("fil_tstart_mjd", hdr["tstart"], fil_path)
    epoch_diff_s = put("inf_epoch_minus_fil_tstart_s",
                       ((imjd - math.floor(tstart)) + (fmjd - (tstart - math.floor(tstart)))) * 86400.0,
                       "inf epoch minus filterbank tstart")
    ok = abs(epoch_diff_s) < dt
    record["checks"]["inf_epoch_matches_fil_tstart_within_one_sample"] = ok
    if not ok:
        raise SystemExit(f"ABORT: .inf epoch and filterbank tstart differ by {epoch_diff_s} s.")

    common = ["prepfold", "-filterbank", "-topo", "-noxwin", "-ncpus", str(args.ncpus),
              "-dm", repr(dm), "-npfact", str(args.npfact), "-ndmfact", str(args.ndmfact)]
    if args.nsub is not None:
        common += ["-nsub", str(args.nsub)]

    # Taylor folds
    for fold, prefix in (("poly", "poly_"), ("kepler_taylor", "kepler_")):
        coeffs = taylor_coeffs(fit, prefix)
        tref = put(f"{fold}_tref_s", fit[f"{prefix}tref_s"][0], fit_path)
        put(f"{fold}_coeffs_at_tref", coeffs, fit_path)
        at_start = put(f"{fold}_coeffs_at_start", reexpand(coeffs, -tref),
                       f"exact re-expansion of {fold}_coeffs_at_tref by -{fold}_tref_s")
        grid = [t_obs * i / 1000.0 for i in range(1001)]
        put(f"{fold}_dropped_phase_at_end_cycles", dropped_phase(at_start, t_obs),
            "derivative orders 3 and above, at t_obs_s")
        put(f"{fold}_dropped_phase_max_abs_cycles",
            max(abs(dropped_phase(at_start, t)) for t in grid),
            "derivative orders 3 and above, max over 1001 points in [0, t_obs_s]")
        f0, f1, f2 = at_start[0], at_start[1], at_start[2]
        cmd = common + ["-f", repr(f0), "-fd", repr(f1), "-fdd", repr(f2)]
        if not args.no_searchfdd:
            cmd.append("-searchfdd")
        cmd += ["-o", f"{args.name}_{fold}", fil_path]
        record["folds"][fold] = shlex.join(cmd)

    # Kepler orbital fold
    f_spin = put("kepler_f_spin_hz", fit["kepler_f_spin_hz"][0], fit_path)
    K = put("kepler_K_hz", fit["kepler_K_hz"][0], fit_path)
    Pb = put("kepler_Pb_s", fit["kepler_Pb_s"][0], fit_path)
    T0 = put("kepler_T0_s", fit["kepler_T0_s"][0], fit_path)
    x = put("orbit_x_lts", K * Pb / (2.0 * math.pi * f_spin),
            "kepler_K_hz * kepler_Pb_s / (2 pi kepler_f_spin_hz)")
    t_asc = put("orbit_T_asc_s", T0 + Pb / 4.0, "kepler_T0_s + kepler_Pb_s / 4")
    to_mjd = put("orbit_To_mjd", imjd + (fmjd + t_asc / 86400.0),
                 "inf epoch + orbit_T_asc_s / 86400")
    cmd = common + ["-bin", "-f", repr(f_spin), "-pb", repr(Pb), "-x", repr(x),
                    "-e", "0", "-w", "0", "-To", f"{to_mjd:.15f}",
                    "-o", f"{args.name}_kepler_orbit", fil_path]
    record["folds"]["kepler_orbit"] = shlex.join(cmd)

    # Write outputs
    os.makedirs(args.outdir, exist_ok=True)
    stem = os.path.join(args.outdir, f"{args.name}_fold_params")
    with open(stem + ".json", "w") as fh:
        json.dump(record, fh, indent=2)
    with open(stem + ".log", "w") as fh:
        fh.write(f"# make_prepfold_inputs.py  {record['created_utc']}\n")
        for k, v in record["inputs"].items():
            fh.write(f"# input {k}: {v}\n")
        for k, v in V.items():
            fh.write(f"{k} = {v['value']!r}    [from: {v['source']}]\n")
        for k, v in record["checks"].items():
            fh.write(f"CHECK {k}: {v}\n")
        for k, v in record["folds"].items():
            fh.write(f"FOLD {k}: {v}\n")
    with open(os.path.join(args.outdir, "prepfold_commands.txt"), "w") as fh:
        for v in record["folds"].values():
            fh.write(v + "\n")

    with open(stem + ".log") as fh:
        print(fh.read())


if __name__ == "__main__":
    main()