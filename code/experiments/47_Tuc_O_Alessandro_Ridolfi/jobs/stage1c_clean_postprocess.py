#!/usr/bin/env python
"""
stage1_clean_postprocess.py

Post-processing for Stage 1 of the 47 Tuc HMM/Viterbi pipeline.

Reads the input and the cleaned filterbank files, computes the per-channel
mean (the "bandpass") using sigpyproc's built-in bandpass() method, produces
a before/after bandpass plot plus zoomed diagnostic plots of known RFI regions.
The zoom plots show BOTH the input and cleaned bandpass, plus an additional
input-only high-resolution zoom to diagnose whether channels between known
RFI bands are genuinely contaminated or clean.

Known RFI sources in MeerKAT L-band (856-1712 MHz), from the SARAO RFI page
(https://skaafrica.atlassian.net/wiki/spaces/ESDKB/pages/305332225/):
    GSM downlink         : 925  - 960  MHz  [hard-zapped]
    Aircraft transponders: 962  - 1213 MHz  [intermittent, narrow -> kadane]
    GPS L2               : 1217 - 1237 MHz  [narrow -> kadane]
    GLONASS L2           : 1242 - 1249 MHz  [narrow -> kadane]
    Inmarsat             : 1526 - 1554 MHz  [hard-zapped]
    GPS L1 / Galileo E1  : 1565 - 1585 MHz  [hard-zapped]
    GLONASS L1           : 1592 - 1610 MHz  [hard-zapped]
    Iridium              : 1616 - 1626 MHz  [hard-zapped]

Inputs
------
    --input-fil      Path to the original (uncleaned) filterbank.
    --cleaned-fil    Path to the cleaned filterbank produced by filtool.
    --stage1-dir     Directory in which to write outputs.
    --filtool-log    Path to the filtool stdout/stderr log (for stats parsing).

Outputs
-------
    <stage1-dir>/bandpass.png                Full-band overview (input + cleaned)
    <stage1-dir>/bandpass_zoom_gsm.png       Zoom: GSM region (880-970 MHz)
    <stage1-dir>/bandpass_zoom_sat.png       Zoom: satellite cluster (1500-1650 MHz)
    <stage1-dir>/bandpass_zoom_iridium.png   Zoom: Iridium + upper edge (1600-1720 MHz)
    <stage1-dir>/bandpass_input_zoom_sat.png Input-only zoom: satellite cluster
                                             to diagnose gaps between RFI bands
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
# Known RFI bands from SARAO L-band documentation.
# Format: (label, flo_MHz, fhi_MHz, hard_zapped)
# hard_zapped=True  -> red  (we hard-zap these in filtool)
# hard_zapped=False -> orange (left to kadane)
# ---------------------------------------------------------------------------
RFI_BANDS = [
    ("GSM downlink",          925.0,  960.0,  True),
    ("Aircraft transponders",  962.0, 1213.0,  False),
    ("GPS L2",                1217.0, 1237.0,  False),
    ("GLONASS L2",            1242.0, 1249.0,  False),
    ("Inmarsat",              1526.0, 1554.0,  True),
    ("GPS L1 / Galileo E1",   1565.0, 1585.0,  True),
    ("GLONASS L1",            1592.0, 1610.0,  True),
    ("Iridium",               1616.0, 1626.0,  True),
]

# Hard-zap ranges actually applied to filtool (derived from RFI_BANDS)
ZAP_BANDS_MHZ = [(flo, fhi) for _, flo, fhi, zapped in RFI_BANDS if zapped]

# Zoom regions for input+cleaned comparison plots
ZOOM_REGIONS = [
    ("GSM region",            880.0,  970.0, "bandpass_zoom_gsm.png"),
    ("Satellite cluster",    1500.0, 1650.0, "bandpass_zoom_sat.png"),
    ("Iridium + upper edge", 1600.0, 1720.0, "bandpass_zoom_iridium.png"),
]

# Input-only zoom regions to inspect gaps between known RFI bands
INPUT_ONLY_ZOOMS = [
    ("Satellite cluster: input only (to diagnose gaps between RFI bands)",
     1500.0, 1650.0, "bandpass_input_zoom_sat.png"),
]


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
    chan_mean  : np.ndarray
    nsamples   : int
    """
    f = FilReader(str(fil_path))
    h = f.header
    nsamples = int(getattr(h, "nsamples", getattr(h, "nspectra", 0)))

    bp = f.bandpass(gulp=gulp)
    bp_arr = np.asarray(bp, dtype=np.float64)
    chan_mean = bp_arr / float(nsamples)

    return channel_freqs(h), chan_mean, nsamples


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def _shade_rfi_bands(ax, flo_view, fhi_view, show_legend=True):
    """Shade known RFI bands within the visible frequency range."""
    seen = {}
    for label, flo, fhi, zapped in RFI_BANDS:
        if fhi < flo_view or flo > fhi_view:
            continue
        color = "red" if zapped else "orange"
        suffix = "zapped" if zapped else "kadane"
        legend_label = f"{label} ({suffix})"
        patch = ax.axvspan(flo, fhi, color=color, alpha=0.15,
                           label=legend_label if legend_label not in seen else "_nolegend_")
        seen[legend_label] = patch


def _dedupe_legend(ax, **kwargs):
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), **kwargs)


def plot_bandpass(in_f, in_mean, cl_f, cl_mean, output_path):
    """Full-band overview: input vs cleaned."""
    fig, ax = plt.subplots(figsize=(12, 4))

    ax.plot(in_f, in_mean, color="black", lw=0.8, label="Input")
    ax.plot(cl_f, cl_mean, color="C0",   lw=0.8, alpha=0.85, label="Cleaned")
    _shade_rfi_bands(ax, float(in_f.min()), float(in_f.max()))

    ax.set_ylabel("Per-channel mean")
    ax.set_xlabel("Frequency (MHz)")
    ax.grid(True, alpha=0.3)
    _dedupe_legend(ax, loc="best", frameon=False, fontsize=7)
    fig.suptitle("Bandpass: input vs cleaned filterbank  "
                 "(red = hard-zapped, orange = kadane only)")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Written: {output_path}")


def plot_zoom(in_f, in_mean, cl_f, cl_mean, flo, fhi, title, output_path):
    """Zoomed input+cleaned plot for a single frequency region."""
    mask_in = (in_f >= flo) & (in_f <= fhi)
    mask_cl = (cl_f >= flo) & (cl_f <= fhi)

    if not mask_in.any():
        print(f"  No channels in zoom region {flo}-{fhi} MHz, skipping.")
        return

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(in_f[mask_in], in_mean[mask_in], color="black", lw=0.8, label="Input")
    ax.plot(cl_f[mask_cl], cl_mean[mask_cl], color="C0",   lw=0.8,
            alpha=0.85, label="Cleaned")
    _shade_rfi_bands(ax, flo, fhi)

    ax.set_xlim(flo, fhi)
    ax.set_ylabel("Per-channel mean")
    ax.set_xlabel("Frequency (MHz)")
    ax.grid(True, alpha=0.3)
    _dedupe_legend(ax, loc="best", frameon=False, fontsize=7)
    fig.suptitle(f"Bandpass zoom: {title}  ({flo:.0f}-{fhi:.0f} MHz)  "
                 "(red = hard-zapped, orange = kadane only)")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Written: {output_path}")


def plot_input_only_zoom(in_f, in_mean, flo, fhi, title, output_path):
    """
    Input-only zoomed plot.

    Shows only the raw (uncleaned) bandpass so we can see whether channels
    between known RFI bands are genuinely elevated above the noise floor,
    independently of what filtool did to them.
    """
    mask = (in_f >= flo) & (in_f <= fhi)

    if not mask.any():
        print(f"  No channels in input-only zoom {flo}-{fhi} MHz, skipping.")
        return

    # Estimate the clean-band baseline from channels outside all RFI bands
    clean_mask = np.ones(len(in_f), dtype=bool)
    for _, blo, bhi, _ in RFI_BANDS:
        clean_mask &= ~((in_f >= blo) & (in_f <= bhi))
    baseline = float(np.median(in_mean[clean_mask])) if clean_mask.any() else float(np.median(in_mean))

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(in_f[mask], in_mean[mask], color="black", lw=0.8, label="Input (raw)")
    ax.axhline(baseline, color="green", lw=1.0, ls="--",
               label=f"Clean-band baseline ({baseline:.2f})")
    _shade_rfi_bands(ax, flo, fhi)

    ax.set_xlim(flo, fhi)
    ax.set_ylabel("Per-channel mean")
    ax.set_xlabel("Frequency (MHz)")
    ax.grid(True, alpha=0.3)
    _dedupe_legend(ax, loc="best", frameon=False, fontsize=7)
    fig.suptitle(f"Input bandpass only: {title}  ({flo:.0f}-{fhi:.0f} MHz)\n"
                 "Green line = clean-band baseline. "
                 "Channels above baseline are RFI-elevated.")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Written: {output_path}")


# ---------------------------------------------------------------------------
# RFI statistics
# ---------------------------------------------------------------------------

def compute_rfi_stats(in_f, in_mean, cl_f, cl_mean):
    with np.errstate(divide="ignore", invalid="ignore"):
        rel_change = np.where(in_mean > 0, (cl_mean - in_mean) / in_mean, 0.0)
    big_drop = rel_change < -0.5

    zapped = np.zeros_like(in_f, dtype=bool)
    for flo, fhi in ZAP_BANDS_MHZ:
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
    parser.add_argument("--input-fil",   required=True)
    parser.add_argument("--cleaned-fil", required=True)
    parser.add_argument("--stage1-dir",  required=True)
    parser.add_argument("--filtool-log", required=True)
    args = parser.parse_args()

    input_fil   = Path(args.input_fil)
    cleaned_fil = Path(args.cleaned_fil)
    stage1_dir  = Path(args.stage1_dir)
    filtool_log = Path(args.filtool_log)
    stage1_dir.mkdir(parents=True, exist_ok=True)

    print(f"Mean bandpass: input   {input_fil}")
    in_f, in_mean, in_n = compute_mean_bandpass(input_fil)
    print(f"Mean bandpass: cleaned {cleaned_fil}")
    cl_f, cl_mean, cl_n = compute_mean_bandpass(cleaned_fil)

    print("\nWriting plots:")
    plot_bandpass(in_f, in_mean, cl_f, cl_mean,
                  stage1_dir / "bandpass.png")

    for title, flo, fhi, fname in ZOOM_REGIONS:
        plot_zoom(in_f, in_mean, cl_f, cl_mean,
                  flo, fhi, title, stage1_dir / fname)

    for title, flo, fhi, fname in INPUT_ONLY_ZOOMS:
        plot_input_only_zoom(in_f, in_mean,
                             flo, fhi, title, stage1_dir / fname)

    stats = compute_rfi_stats(in_f, in_mean, cl_f, cl_mean)
    stats["nsamples_input"]   = int(in_n)
    stats["nsamples_cleaned"] = int(cl_n)
    stats["rfi_flags"] = {
        "algorithms":    ["zdot", "kadaneF 8 4", "kadaneT 8 4"],
        "hard_zaps_mhz": [list(z) for z in ZAP_BANDS_MHZ],
        "fill_patch":    "rand",
        "sarao_rfi_reference": (
            "https://skaafrica.atlassian.net/wiki/spaces/ESDKB/"
            "pages/305332225/Radio+Frequency+Interference+RFI"
        ),
    }
    stats["filtool_log_excerpt"] = parse_filtool_log(filtool_log)
    stats["postprocess_datetime_utc"] = datetime.datetime.now(
        datetime.timezone.utc).isoformat()
    stats["input_fil"]   = str(input_fil)
    stats["cleaned_fil"] = str(cleaned_fil)

    stats_yaml = stage1_dir / "rfi_stats.yaml"
    print(f"\nWriting stats: {stats_yaml}")
    with open(stats_yaml, "w") as fh:
        fh.write("# Stage 1 RFI cleaning statistics.\n"
                 "# Auto-generated by stage1_clean_postprocess.py.\n\n")
        yaml.safe_dump(stats, fh, sort_keys=False, default_flow_style=False)

    print("Stage 1 post-processing complete.")


if __name__ == "__main__":
    main()