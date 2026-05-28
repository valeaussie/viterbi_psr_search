#!/usr/bin/env python
"""
stage0_inspect.py

Stage 0 of the 47 Tuc HMM/Viterbi pulsar search pipeline.

Reads a filterbank file, extracts metadata, computes derived quantities
(Nyquist frequency, Chandler ceiling for 47 Tuc O as a worked example),
queries the ATNF pulsar catalogue for all known 47 Tucanae pulsars
(prefix J0024-7204), computes per-pulsar Chandler ceilings where the
orbital parameters are available, and writes:

    <exp-dir>/config/observation.yaml
    <exp-dir>/config/known_47tuc_pulsars.yaml
    <exp-dir>/stage0_inspect/inspect.log
    <exp-dir>/stage0_inspect/header_raw.txt
    <exp-dir>/provenance/stage0_runinfo.txt

Usage
-----
    python stage0_inspect.py --fil <path-to-fil> --exp-dir <experiment-dir>
    python stage0_inspect.py --fil <path-to-fil> --exp-dir <experiment-dir> --offline

The --offline flag skips the ATNF query (useful on compute nodes
without network access).

Dependencies
------------
    sigpyproc, psrqpy, numpy, pyyaml
"""

import argparse
import datetime
import math
import os
import socket
import subprocess
import sys
from pathlib import Path

import numpy as np
import yaml

# ---------------------------------------------------------------------------
# Constants and reference values
# ---------------------------------------------------------------------------

C_LIGHT = 299792458.0  # m/s, speed of light (exact, SI)

# 47 Tuc O reference values for worked-example Chandler ceiling
TUC_O_F0_HZ = 378.3087883600985037
TUC_O_A1_LTS = 0.045153297
TUC_O_PB_S = 0.135974305889033 * 86400.0  # days -> seconds


# ---------------------------------------------------------------------------
# Filterbank header reading
# ---------------------------------------------------------------------------

def read_filterbank_header(fil_path):
    """
    Read filterbank header using sigpyproc.FilReader.

    Returns a dict of metadata used by downstream stages.
    """
    from sigpyproc import FilReader

    f = FilReader(str(fil_path))
    h = f.header

    # sigpyproc field names depend on version; we try the common ones.
    def get(*names, default=None):
        for n in names:
            if hasattr(h, n):
                return getattr(h, n)
        return default

    meta = {
        "filepath": str(Path(fil_path).resolve()),
        "source_name": get("source_name", "source"),
        "telescope": get("telescope", "telescope_id"),
        "mjd_start": float(get("tstart", "mjd_start")),
        "ra_j2000": str(get("ra", "src_raj")),
        "dec_j2000": str(get("dec", "src_dej")),
        "tsamp_s": float(get("tsamp")),
        "nsamples": int(get("nsamples", "nspectra")),
        "tobs_s": float(get("tobs")),
        "nchan": int(get("nchans", "nchan")),
        "fch1_mhz": float(get("fch1")),
        "foff_mhz": float(get("foff")),
        "nbits": int(get("nbits", default=0)),
        "npol": int(get("nifs", "npol", default=0)),
    }

    # Derived band edges
    fch1 = meta["fch1_mhz"]
    foff = meta["foff_mhz"]
    nchan = meta["nchan"]
    band_lo = fch1 + foff * (nchan - 1) if foff < 0 else fch1
    band_hi = fch1 if foff < 0 else fch1 + foff * (nchan - 1)
    meta["flo_mhz"] = float(min(band_lo, band_hi))
    meta["fhi_mhz"] = float(max(band_lo, band_hi))
    meta["bandwidth_mhz"] = float(abs(foff) * nchan)
    meta["fcentre_mhz"] = float(0.5 * (meta["flo_mhz"] + meta["fhi_mhz"]))

    return meta, h


def header_to_text(h):
    """Dump every attribute of the sigpyproc header for the raw log."""
    lines = []
    for attr in sorted(dir(h)):
        if attr.startswith("_"):
            continue
        try:
            value = getattr(h, attr)
        except Exception as exc:
            value = f"<error: {exc}>"
        if callable(value):
            continue
        lines.append(f"{attr} = {value}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------------------------

def sanity_checks(meta):
    """
    Check internal consistency of header values.

    Returns a dict with results, suitable for embedding in YAML.
    """
    checks = {}

    # tsamp * nsamples ~ tobs
    computed_tobs = meta["tsamp_s"] * meta["nsamples"]
    residual = abs(computed_tobs - meta["tobs_s"])
    checks["tsamp_nsamples_consistent"] = bool(residual < meta["tsamp_s"])
    checks["tsamp_nsamples_residual_s"] = float(residual)
    checks["computed_tobs_s"] = float(computed_tobs)

    # File-on-disk consistency
    if Path(meta["filepath"]).exists():
        checks["file_exists"] = True
        checks["file_size_bytes"] = int(Path(meta["filepath"]).stat().st_size)
    else:
        checks["file_exists"] = False
        checks["file_size_bytes"] = None

    return checks


# ---------------------------------------------------------------------------
# Derived quantities
# ---------------------------------------------------------------------------

def chandler_tcoh_max(f_p_hz, p_b_s, a1_lts):
    """
    Chandler (2003) ceiling on Tcoh for a nearly circular binary.

        T_coh_max = sqrt( P_b / (2 pi f_p beta) )

    with beta = (V_1 sin i)/c = 2 pi (a_1 sin i) / (c P_b).

    Parameters
    ----------
    f_p_hz : float
        Pulse frequency (Hz).
    p_b_s : float
        Orbital period (seconds).
    a1_lts : float
        Projected semi-major axis (light-seconds).

    Returns
    -------
    float
        T_coh_max in seconds.
    """
    # a1_lts is already in light-seconds, so V_1 sin i / c = 2 pi a1_lts / P_b
    beta = 2.0 * math.pi * a1_lts / p_b_s
    return math.sqrt(p_b_s / (2.0 * math.pi * f_p_hz * beta))


def derived_quantities(meta):
    """Compute derived quantities used in the pipeline planning."""
    nyq = 0.5 / meta["tsamp_s"]
    return {
        "nyquist_hz": float(nyq),
        "recommended_search_flo_hz": 10.0,
        "recommended_search_fhi_hz": 1000.0,
        "tuc_o_chandler_tcoh_max_s": float(
            chandler_tcoh_max(TUC_O_F0_HZ, TUC_O_PB_S, TUC_O_A1_LTS)
        ),
    }


# ---------------------------------------------------------------------------
# ATNF query for 47 Tuc pulsars
# ---------------------------------------------------------------------------

def query_atnf_47tuc():
    """
    Query the ATNF pulsar catalogue for all pulsars matching J0024-7204*.

    Returns
    -------
    list of dict
        One entry per pulsar with available parameters.
    str
        ATNF catalogue version string.
    """
    from psrqpy import QueryATNF

    params = [
        "JNAME", "PSRJ", "RAJ", "DECJ",
        "P0", "F0", "DM",
        "PB", "A1", "ECC", "OM", "T0",
        "BINARY",
    ]

    q = QueryATNF(params=params)
    df = q.dataframe

    # Filter to 47 Tuc by name prefix.
    name_col = "PSRJ" if "PSRJ" in df.columns else "JNAME"
    mask = df[name_col].astype(str).str.startswith("J0024-7204")
    sub = df[mask].copy()

    pulsars = []
    for _, row in sub.iterrows():
        entry = {
            "name": str(row.get(name_col, "")),
            "raj": _to_str(row.get("RAJ")),
            "decj": _to_str(row.get("DECJ")),
            "P0_s": _to_float(row.get("P0")),
            "F0_hz": _to_float(row.get("F0")),
            "DM_pccm3": _to_float(row.get("DM")),
            "PB_days": _to_float(row.get("PB")),
            "A1_lts": _to_float(row.get("A1")),
            "ECC": _to_float(row.get("ECC")),
            "OM_deg": _to_float(row.get("OM")),
            "T0_mjd": _to_float(row.get("T0")),
            "BINARY_type": _to_str(row.get("BINARY")) if row.get("BINARY") else None,
        }

        # Compute Chandler ceiling where possible.
        if (entry["F0_hz"] is not None and
                entry["PB_days"] is not None and
                entry["A1_lts"] is not None):
            pb_s = entry["PB_days"] * 86400.0
            entry["Tcoh_max_s"] = float(
                chandler_tcoh_max(entry["F0_hz"], pb_s, entry["A1_lts"])
            )
            entry["is_binary"] = True
        elif entry["PB_days"] is None and entry["F0_hz"] is not None:
            entry["Tcoh_max_s"] = None  # isolated -> no Doppler ceiling
            entry["is_binary"] = False
        else:
            entry["Tcoh_max_s"] = None  # incomplete entry
            entry["is_binary"] = None

        pulsars.append(entry)

    pulsars.sort(key=lambda x: x["name"])

    try:
        version = str(q.get_catalogue_version())
    except Exception:
        version = "unknown"

    return pulsars, version


def _to_float(v):
    if v is None:
        return None
    try:
        if isinstance(v, float) and math.isnan(v):
            return None
    except Exception:
        pass
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _to_str(v):
    if v is None:
        return None
    try:
        if isinstance(v, float) and math.isnan(v):
            return None
    except Exception:
        pass
    return str(v)


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

def git_hash(path):
    """Return the git hash of the working directory, or 'unknown'."""
    try:
        out = subprocess.check_output(
            ["git", "-C", str(path), "rev-parse", "HEAD"],
            stderr=subprocess.DEVNULL,
        )
        return out.decode().strip()
    except Exception:
        return "unknown"


def make_provenance(script_path):
    return {
        "script": str(Path(script_path).resolve()),
        "git_hash": git_hash(Path(script_path).parent),
        "hostname": socket.gethostname(),
        "user": os.environ.get("USER", "unknown"),
        "datetime_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "python": sys.version.replace("\n", " "),
        "argv": sys.argv,
    }


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def ensure_dirs(exp_dir):
    exp_dir = Path(exp_dir)
    (exp_dir / "config").mkdir(parents=True, exist_ok=True)
    (exp_dir / "stage0_inspect").mkdir(parents=True, exist_ok=True)
    (exp_dir / "provenance").mkdir(parents=True, exist_ok=True)
    return exp_dir


def write_yaml(path, data, header_comment=None):
    with open(path, "w") as fh:
        if header_comment:
            for line in header_comment.splitlines():
                fh.write(f"# {line}\n")
            fh.write("\n")
        yaml.safe_dump(data, fh, sort_keys=False, default_flow_style=False)


def fmt_pulsar_line(p):
    """One-line summary for a pulsar entry."""
    f0 = f"{p['F0_hz']:.4f}" if p.get("F0_hz") is not None else "  -   "
    p0 = f"{p['P0_s']*1e3:.3f}" if p.get("P0_s") is not None else "  -  "
    dm = f"{p['DM_pccm3']:.3f}" if p.get("DM_pccm3") is not None else "  -  "
    pb = f"{p['PB_days']:.4f}" if p.get("PB_days") is not None else "   -   "
    a1 = f"{p['A1_lts']:.5f}" if p.get("A1_lts") is not None else "   -   "
    tch = (
        f"{p['Tcoh_max_s']:.1f}"
        if p.get("Tcoh_max_s") is not None
        else ("inf" if p.get("is_binary") is False else "  -  ")
    )
    return (f"  {p['name']:<14}  P0={p0:>7} ms  DM={dm:>7}  "
            f"PB={pb:>9} d  A1={a1:>9} lt-s  Tcoh_max={tch:>8} s")


def write_log(path, meta, sanity, derived, pulsars, atnf_version, offline):
    lines = []
    lines.append("=" * 78)
    lines.append("Stage 0 inspection report")
    lines.append(datetime.datetime.now(datetime.timezone.utc).isoformat())
    lines.append("=" * 78)

    lines.append("")
    lines.append("Observation")
    lines.append("-----------")
    for k in ["filepath", "source_name", "telescope", "mjd_start",
              "ra_j2000", "dec_j2000"]:
        lines.append(f"  {k:>20s} : {meta[k]}")

    lines.append("")
    lines.append("Sampling and band")
    lines.append("-----------------")
    for k in ["tsamp_s", "nsamples", "tobs_s", "nchan",
              "flo_mhz", "fhi_mhz", "bandwidth_mhz", "fcentre_mhz",
              "nbits", "npol"]:
        lines.append(f"  {k:>20s} : {meta[k]}")

    lines.append("")
    lines.append("Sanity checks")
    lines.append("-------------")
    for k, v in sanity.items():
        lines.append(f"  {k:>30s} : {v}")

    lines.append("")
    lines.append("Derived quantities")
    lines.append("------------------")
    for k, v in derived.items():
        lines.append(f"  {k:>30s} : {v}")

    lines.append("")
    lines.append(f"Known 47 Tuc pulsars  (ATNF catalogue version: {atnf_version})")
    lines.append("-" * 78)
    if offline:
        lines.append("  [skipped: --offline mode]")
    elif not pulsars:
        lines.append("  [no entries returned]")
    else:
        lines.append(f"  {len(pulsars)} pulsars matched J0024-7204*:")
        for p in pulsars:
            lines.append(fmt_pulsar_line(p))

    lines.append("")
    lines.append("=" * 78)

    text = "\n".join(lines)
    with open(path, "w") as fh:
        fh.write(text + "\n")
    return text


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Stage 0: inspect filterbank metadata for the 47 Tuc HMM/Viterbi pipeline."
    )
    parser.add_argument("--fil", required=True,
                        help="Path to the input filterbank (.fil) file.")
    parser.add_argument("--exp-dir", required=True,
                        help="Experiment directory (will be created if missing).")
    parser.add_argument("--offline", action="store_true",
                        help="Skip the ATNF query (no network).")
    args = parser.parse_args()

    fil_path = Path(args.fil)
    if not fil_path.exists():
        sys.exit(f"ERROR: input file not found: {fil_path}")

    exp_dir = ensure_dirs(args.exp_dir)

    # 1. Header
    meta, h = read_filterbank_header(fil_path)

    # 2. Sanity
    sanity = sanity_checks(meta)

    # 3. Derived
    derived = derived_quantities(meta)

    # 4. ATNF query
    if args.offline:
        pulsars, atnf_version = [], "skipped"
    else:
        try:
            pulsars, atnf_version = query_atnf_47tuc()
        except Exception as exc:
            print(f"WARNING: ATNF query failed: {exc}", file=sys.stderr)
            pulsars, atnf_version = [], f"failed: {exc}"

    # 5. Write outputs
    obs_yaml = exp_dir / "config" / "observation.yaml"
    write_yaml(
        obs_yaml,
        {
            "observation": {
                "filepath": meta["filepath"],
                "source_name": meta["source_name"],
                "telescope": meta["telescope"],
                "mjd_start": meta["mjd_start"],
                "ra_j2000": meta["ra_j2000"],
                "dec_j2000": meta["dec_j2000"],
            },
            "sampling": {
                "tsamp_s": meta["tsamp_s"],
                "nsamples": meta["nsamples"],
                "tobs_s": meta["tobs_s"],
                "nchan": meta["nchan"],
                "fcentre_mhz": meta["fcentre_mhz"],
                "bandwidth_mhz": meta["bandwidth_mhz"],
                "flo_mhz": meta["flo_mhz"],
                "fhi_mhz": meta["fhi_mhz"],
                "nbits": meta["nbits"],
                "npol": meta["npol"],
            },
            "derived": derived,
            "sanity_checks": sanity,
        },
        header_comment=(
            "Observation metadata for the 47 Tuc HMM/Viterbi pipeline.\n"
            "Auto-generated by stage0_inspect.py. Do not edit by hand."
        ),
    )

    pulsars_yaml = exp_dir / "config" / "known_47tuc_pulsars.yaml"
    write_yaml(
        pulsars_yaml,
        {
            "atnf_query": {
                "catalogue_version": atnf_version,
                "query_datetime_utc": datetime.datetime.now(
                    datetime.timezone.utc
                ).isoformat(),
                "name_prefix": "J0024-7204",
                "offline": bool(args.offline),
            },
            "pulsars": pulsars,
        },
        header_comment=(
            "Known 47 Tuc pulsars (ATNF catalogue, J0024-7204* prefix).\n"
            "Per-pulsar Chandler ceiling T_coh_max computed by stage0_inspect.py.\n"
            "Isolated pulsars: T_coh_max = None (unbounded).\n"
            "Pulsars with incomplete orbital info: T_coh_max = None and is_binary = None."
        ),
    )

    header_dump = exp_dir / "stage0_inspect" / "header_raw.txt"
    with open(header_dump, "w") as fh:
        fh.write(header_to_text(h))

    log_path = exp_dir / "stage0_inspect" / "inspect.log"
    text = write_log(log_path, meta, sanity, derived, pulsars, atnf_version,
                     args.offline)

    prov = make_provenance(__file__)
    prov_path = exp_dir / "provenance" / "stage0_runinfo.txt"
    with open(prov_path, "w") as fh:
        yaml.safe_dump(prov, fh, sort_keys=False, default_flow_style=False)

    # 6. Echo summary to stdout
    print(text)
    print(f"\nWrote:\n  {obs_yaml}\n  {pulsars_yaml}\n"
          f"  {header_dump}\n  {log_path}\n  {prov_path}")


if __name__ == "__main__":
    main()