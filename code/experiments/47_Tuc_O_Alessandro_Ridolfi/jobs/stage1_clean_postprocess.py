#!/usr/bin/env python
"""
stage1_clean_postprocess.py

Post-processing for Stage 1 of the 47 Tuc HMM/Viterbi pipeline.

Reads the input and the cleaned filterbank files, computes the per-channel
mean (the "bandpass") using sigpyproc's built-in bandpass() method, optionally
estimates a per-channel median from a subsample of blocks, produces a
before/after bandpass plot, and writes RFI statistics to YAML.

Inputs
------
    --input-fil      Path to the original (uncleaned) filterbank.
    --cleaned-fil    Path to the cleaned filterbank produced by filtool.
    --stage1-dir     Directory in which to write outputs.
    --filtool-log    Path to the filtool stdout/stderr log (for stats parsing).
    --no-median      Skip the (slower) per-channel median estimate.

Outputs
-------
    <stage1-dir>/bandpass.png
    <stage1-dir>/rfi_stats.yaml
"""

import argparse
import datetime
from pathlib import Path

import numpy as np
import yaml

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sigpyproc import FilReader


# ---------------------------------------------------------------------------
# Bandpass computation
# ---------------------------------------------------------------------------

def channel_freqs(h):
    """Channel centre frequencies (MHz) from a sigpyproc header."""
    nchan = int(getattr(h, "nchans", getattr(h, "nchan", 0)))
    fch1 = float(h.fch1)
    foff = float(h.foff)
    return fch1 + foff * np.arange(nchan)


def compute_mean_bandpass(fil_path, gulp=4096):
    """
    Per-channel mean intensity over the whole observation.

    Uses sigpyproc's bandpass() (sum across time), then divides by the number
    of samples to get the mean.

    Returns
    -------
    freqs_mhz : np.ndarray
    chan_mean : np.ndarray
    nsamples : int
    """
    f = FilReader(str(fil_path))
    h = f.header
    nsamples = int(getattr(h, "nsamples", getattr(h, "nspectra", 0)))

    bp = f.bandpass(gulp=gulp)          # TimeSeries of length nchan (sum over time)
    bp_arr = np.asarray(bp, dtype=np.float64)
    chan_mean = bp_arr / float(nsamples)

    return channel_freqs(h), chan_mean, nsamples


def compute_median_bandpass(fil_path, n_blocks=20, gulp=65536, seed=0):
    """
    Per-channel median estimate from a random subsample of blocks.

    Reads n_blocks evenly spaced blocks of gulp samples each and concatenates
    them for the median. This is an estimate, not the exact full-file median,
    but is accurate to a fraction of a per cent for diagnostic purposes.

    Returns
    -------
    chan_median : np.ndarray
    """
    f = FilReader(str(fil_path))
    h = f.header
    nchan = int(getattr(h, "nchans", getattr(h, "nchan", 0)))
    nsamples = int(getattr(h, "nsamples", getattr(h, "nspectra", 0)))

    max_start = max(nsamples - gulp, 1)
    starts = np.linspace(0, max_start, n_blocks, dtype=int)

    chunks = []
    for s in starts:
        nread = min(gulp, nsamples - int(s))
        if nread <= 0:
            continue
        block = f.readBlock(int(s), int(nread))
        arr = np.asarray(block, dtype=np.float32)   # (nchan, nread) expected
        if arr.shape[0] != nchan and arr.shape[-1] == nchan:
            arr = arr.T
        chunks.append(arr)

    if not chunks:
        return np.full(nchan, np.nan)

    alldata = np.concatenate(chunks, axis=1)        # (nchan, total_samples)
    return np.median(alldata, axis=1)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_bandpass(in_f, in_mean, in_med, cl_f, cl_mean, cl_med,
                  output_path, zap_bands_mhz, have_median):
    nrows = 2 if have_median else 1
    fig, axes = plt.subplots(nrows, 1, figsize=(10, 3.5 * nrows), sharex=True)
    if nrows == 1:
        axes = [axes]

    ax = axes[0]
    ax.plot(in_f, in_mean, color="black", lw=0.8, label="Input")
    ax.plot(cl_f, cl_mean, color="C0", lw=0.8, alpha=0.85, label="Cleaned")
    for flo, fhi in zap_bands_mhz:
        ax.axvspan(flo, fhi, color="red", alpha=0.12)
    ax.set_ylabel("Per-channel mean")
    ax.legend(loc="best", frameon=False)
    ax.grid(True, alpha=0.3)

    if have_median:
        ax = axes[1]
        ax.plot(in_f, in_med, color="black", lw=0.8, label="Input")
        ax.plot(cl_f, cl_med, color="C0", lw=0.8, alpha=0.85, label="Cleaned")
        for flo, fhi in zap_bands_mhz:
            ax.axvspan(flo, fhi, color="red", alpha=0.12)
        ax.set_ylabel("Per-channel median")
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel("Frequency (MHz)")
    fig.suptitle("Bandpass: input vs cleaned filterbank "
                 "(red bands = hard-zapped)")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# RFI statistics
# ---------------------------------------------------------------------------

def compute_rfi_stats(in_f, in_mean, cl_f, cl_mean, zap_bands_mhz):
    with np.errstate(divide="ignore", invalid="ignore"):
        rel_change = np.where(in_mean > 0, (cl_mean - in_mean) / in_mean, 0.0)
    big_drop = rel_change < -0.5

    zapped = np.zeros_like(in_f, dtype=bool)
    for flo, fhi in zap_bands_mhz:
        zapped |= (in_f >= flo) & (in_f <= fhi)

    return {
        "channels_total": int(in_f.size),
        "channels_strongly_affected_by_cleaning": int(big_drop.sum()),
        "fraction_strongly_affected_by_cleaning": float(big_drop.mean()),
        "fraction_band_hard_zapped": float(zapped.mean()),
        "mean_intensity_input_global": float(in_mean.mean()),
        "mean_intensity_cleaned_global": float(cl_mean.mean()),
    }


# ---------------------------------------------------------------------------
# filtool log parsing (best-effort)
# ---------------------------------------------------------------------------

def parse_filtool_log(log_path):
    if not log_path.exists():
        return {"log_present": False}
    keywords = ["flag", "zap", "mask", "rfi", "kadane", "segment", "patch", "fill"]
    matches = []
    with open(log_path, "r", errors="replace") as fh:
        for line in fh:
            if any(k in line.lower() for k in keywords):
                matches.append(line.rstrip())
    return {"log_present": True, "log_path": str(log_path),
            "lines_matching_rfi_keywords": matches[-200:]}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-fil", required=True)
    parser.add_argument("--cleaned-fil", required=True)
    parser.add_argument("--stage1-dir", required=True)
    parser.add_argument("--filtool-log", required=True)
    parser.add_argument("--no-median", action="store_true",
                        help="Skip the per-channel median estimate.")
    args = parser.parse_args()

    input_fil = Path(args.input_fil)
    cleaned_fil = Path(args.cleaned_fil)
    stage1_dir = Path(args.stage1_dir)
    filtool_log = Path(args.filtool_log)
    stage1_dir.mkdir(parents=True, exist_ok=True)

    zap_bands_mhz = [(925.0, 960.0), (1525.0, 1612.0), (1675.0, 1720.0)]

    print(f"Mean bandpass: input  {input_fil}")
    in_f, in_mean, in_n = compute_mean_bandpass(input_fil)
    print(f"Mean bandpass: cleaned {cleaned_fil}")
    cl_f, cl_mean, cl_n = compute_mean_bandpass(cleaned_fil)

    have_median = not args.no_median
    if have_median:
        print("Median bandpass: input (subsample)")
        in_med = compute_median_bandpass(input_fil)
        print("Median bandpass: cleaned (subsample)")
        cl_med = compute_median_bandpass(cleaned_fil)
    else:
        in_med = cl_med = None

    png = stage1_dir / "bandpass.png"
    print(f"Writing plot: {png}")
    plot_bandpass(in_f, in_mean, in_med, cl_f, cl_mean, cl_med,
                  png, zap_bands_mhz, have_median)

    stats = compute_rfi_stats(in_f, in_mean, cl_f, cl_mean, zap_bands_mhz)
    stats["nsamples_input"] = int(in_n)
    stats["nsamples_cleaned"] = int(cl_n)
    stats["rfi_flags"] = {
        "algorithms": ["zdot", "kadaneF 8 4", "kadaneT 8 4"],
        "hard_zaps_mhz": [list(z) for z in zap_bands_mhz],
        "fill_patch": "rand",
    }
    stats["filtool_log_excerpt"] = parse_filtool_log(filtool_log)
    stats["postprocess_datetime_utc"] = datetime.datetime.now(
        datetime.timezone.utc).isoformat()
    stats["input_fil"] = str(input_fil)
    stats["cleaned_fil"] = str(cleaned_fil)

    stats_yaml = stage1_dir / "rfi_stats.yaml"
    print(f"Writing stats: {stats_yaml}")
    with open(stats_yaml, "w") as fh:
        fh.write("# Stage 1 RFI cleaning statistics.\n"
                 "# Auto-generated by stage1_clean_postprocess.py.\n\n")
        yaml.safe_dump(stats, fh, sort_keys=False, default_flow_style=False)

    print("Stage 1 post-processing complete.")


if __name__ == "__main__":
    main()