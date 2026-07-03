#!/usr/bin/env python3
"""
inspect_dynamic_spectrum.py

Plots dynamic spectra (time vs frequency, colour = intensity) for a short
chunk of a filterbank file, focused on a specified frequency range.
Also computes per-channel statistics (mean, std, max) for each chunk and
compares RFI band channels against the clean-band baseline.

Run this on the RAW (uncleaned) filterbank to see whether RFI appears as
bright streaks at known satellite frequencies.

Three chunks are plotted: start, middle, and one-third of the way through
the observation, to catch both persistent and intermittent RFI.

Usage
-----
    python inspect_dynamic_spectrum.py \
        --fil   path/to/raw.fil \
        --flo   1500 \
        --fhi   1650 \
        --nsamples 10000 \
        --outdir path/to/stage1_clean/

Outputs
-------
    <outdir>/dynspec_start.png
    <outdir>/dynspec_middle.png
    <outdir>/dynspec_onethird.png
    <outdir>/dynspec_stats.txt     Per-channel statistics for each chunk
    <outdir>/dynspec_rfi_ratio.png Per-channel mean/baseline ratio plot
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sigpyproc import FilReader


# Known RFI bands for reference shading (same as stage1_clean_postprocess.py)
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


def channel_freqs(h):
    nchan = int(getattr(h, "nchans", getattr(h, "nchan", 0)))
    fch1  = float(h.fch1)
    foff  = float(h.foff)
    return fch1 + foff * np.arange(nchan)


def read_chunk(fil_path, start_sample, nsamples):
    f = FilReader(str(fil_path))
    block = f.readBlock(int(start_sample), int(nsamples))
    arr = np.asarray(block, dtype=np.float32)
    h = f.header
    nchan = int(getattr(h, "nchans", getattr(h, "nchan", 0)))
    if arr.shape[0] != nchan and arr.shape[-1] == nchan:
        arr = arr.T
    return arr   # shape (nchan, nsamples)


def clean_band_mask(freqs):
    """Boolean mask: True for channels outside all known RFI bands."""
    mask = np.ones(len(freqs), dtype=bool)
    for _, blo, bhi, _ in RFI_BANDS:
        mask &= ~((freqs >= blo) & (freqs <= bhi))
    return mask


def compute_channel_stats(freqs, data, flo, fhi):
    """
    For each channel in flo-fhi, compute mean, std, max.
    Also compute the clean-band baseline (median of per-channel means
    outside all RFI bands).

    Returns dict with arrays indexed to the selected frequency range.
    """
    mask_range = (freqs >= flo) & (freqs <= fhi)
    f_sel  = freqs[mask_range]
    d_sel  = data[mask_range, :]   # (n_sel_chans, nsamples)

    chan_mean = d_sel.mean(axis=1)
    chan_std  = d_sel.std(axis=1)
    chan_max  = d_sel.max(axis=1)

    # Baseline from clean channels across the FULL band
    cb_mask   = clean_band_mask(freqs)
    baseline  = float(np.median(data[cb_mask, :].mean(axis=1)))

    return {
        "freqs":    f_sel,
        "mean":     chan_mean,
        "std":      chan_std,
        "max":      chan_max,
        "baseline": baseline,
        "ratio":    chan_mean / baseline,   # >1 means RFI-elevated
    }


def plot_dynspec(freqs, data, flo, fhi, start_sample, tsamp,
                 label, output_path):
    mask = (freqs >= flo) & (freqs <= fhi)
    if not mask.any():
        print(f"  No channels in {flo}-{fhi} MHz, skipping {label}.")
        return

    f_sel    = freqs[mask]
    d_sel    = data[mask, :]
    nsamples = d_sel.shape[1]
    times    = np.arange(nsamples) * tsamp

    vmin = np.percentile(d_sel, 1)
    vmax = np.percentile(d_sel, 99)

    fig, ax = plt.subplots(figsize=(12, 5))
    extent = [times[0], times[-1], f_sel.min(), f_sel.max()]
    im = ax.imshow(
        d_sel, aspect="auto", origin="lower", extent=extent,
        vmin=vmin, vmax=vmax, cmap="plasma", interpolation="nearest",
    )
    plt.colorbar(im, ax=ax, label="Intensity (raw counts)")

    for band_label, blo, bhi, zapped in RFI_BANDS:
        if bhi < flo or blo > fhi:
            continue
        color = "cyan" if zapped else "yellow"
        ax.axhspan(max(blo, flo), min(bhi, fhi), color=color, alpha=0.25,
                   label=f"{band_label} ({'zapped' if zapped else 'kadane'})")

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(),
              loc="upper right", frameon=True, fontsize=7, framealpha=0.7)

    t_start_s = start_sample * tsamp
    ax.set_xlabel("Time within chunk (s)")
    ax.set_ylabel("Frequency (MHz)")
    ax.set_title(
        f"Dynamic spectrum: {label}\n"
        f"Chunk start = {t_start_s:.1f} s  "
        f"({nsamples} samples, {nsamples * tsamp:.2f} s)  "
        f"Freq: {flo}-{fhi} MHz\n"
        f"Cyan = hard-zapped bands, Yellow = kadane-only bands"
    )
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Written: {output_path}")


def plot_rfi_ratio(all_stats, flo, fhi, output_path):
    """
    Plot per-channel mean/baseline ratio for all chunks on one figure.
    Ratio > 1 means the channel mean is elevated above the clean baseline.
    A horizontal line at 1.0 is the baseline.
    A horizontal line at 1.05 (5% above baseline) is a conservative
    threshold: anything above this is worth investigating.
    """
    fig, ax = plt.subplots(figsize=(12, 5))

    colors = {"start": "C0", "onethird": "C1", "middle": "C2"}
    for label, stats in all_stats.items():
        ax.plot(stats["freqs"], stats["ratio"],
                lw=0.8, alpha=0.8, color=colors.get(label, "black"),
                label=f"{label} (baseline={stats['baseline']:.2f})")

    ax.axhline(1.0,  color="black", lw=1.0, ls="--", label="Baseline (ratio=1)")
    ax.axhline(1.05, color="red",   lw=0.8, ls=":",
               label="5% above baseline")

    for band_label, blo, bhi, zapped in RFI_BANDS:
        if bhi < flo or blo > fhi:
            continue
        color = "red" if zapped else "orange"
        ax.axvspan(max(blo, flo), min(bhi, fhi), color=color, alpha=0.12,
                   label=f"{band_label} ({'zapped' if zapped else 'kadane'})")

    handles, labels_leg = ax.get_legend_handles_labels()
    by_label = dict(zip(labels_leg, handles))
    ax.legend(by_label.values(), by_label.keys(),
              loc="upper right", frameon=True, fontsize=7, framealpha=0.8)

    ax.set_xlim(flo, fhi)
    ax.set_xlabel("Frequency (MHz)")
    ax.set_ylabel("Per-channel mean / clean-band baseline")
    ax.set_title(
        "RFI elevation ratio: per-channel mean divided by clean-band baseline\n"
        "Ratio > 1.05 (red dotted line) indicates RFI-elevated channels.\n"
        "Red shading = hard-zapped bands, Orange = kadane-only bands."
    )
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Written: {output_path}")


def write_stats_txt(all_stats, flo, fhi, output_path):
    """
    Write a text table of per-channel statistics for channels inside
    known RFI bands, showing whether they are elevated above baseline.
    """
    lines = []
    lines.append("Per-channel RFI elevation statistics")
    lines.append(f"Frequency range: {flo}-{fhi} MHz")
    lines.append("Ratio = per-channel mean / clean-band baseline")
    lines.append("Channels with ratio > 1.05 are RFI-elevated.")
    lines.append("")

    for band_label, blo, bhi, zapped in RFI_BANDS:
        if bhi < flo or blo > fhi:
            continue
        lines.append(f"{'='*70}")
        lines.append(f"Band: {band_label}  ({blo}-{bhi} MHz)  "
                     f"{'[HARD-ZAPPED]' if zapped else '[kadane only]'}")
        lines.append(f"{'Chunk':<12} {'Baseline':>10} {'Mean in band':>14} "
                     f"{'Ratio':>8} {'Std in band':>13} {'Max in band':>13}")
        lines.append("-" * 70)
        for label, stats in all_stats.items():
            band_mask = (stats["freqs"] >= blo) & (stats["freqs"] <= bhi)
            if not band_mask.any():
                continue
            mean_band = float(stats["mean"][band_mask].mean())
            std_band  = float(stats["std"][band_mask].mean())
            max_band  = float(stats["max"][band_mask].max())
            ratio     = mean_band / stats["baseline"]
            lines.append(
                f"{label:<12} {stats['baseline']:>10.3f} {mean_band:>14.3f} "
                f"{ratio:>8.4f} {std_band:>13.3f} {max_band:>13.1f}"
            )
        lines.append("")

    with open(output_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"  Written: {output_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fil",      required=True)
    parser.add_argument("--flo",      type=float, default=1500.0)
    parser.add_argument("--fhi",      type=float, default=1650.0)
    parser.add_argument("--nsamples", type=int,   default=10000)
    parser.add_argument("--outdir",   required=True)
    args = parser.parse_args()

    fil_path = Path(args.fil)
    outdir   = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    f = FilReader(str(fil_path))
    h = f.header
    nsamples_total = int(getattr(h, "nsamples", getattr(h, "nspectra", 0)))
    tsamp          = float(h.tsamp)
    freqs          = channel_freqs(h)

    print(f"File:           {fil_path}")
    print(f"Total samples:  {nsamples_total}")
    print(f"Tsamp:          {tsamp:.6e} s")
    print(f"Duration:       {nsamples_total * tsamp / 3600:.2f} hours")
    print(f"Chunk size:     {args.nsamples} samples = "
          f"{args.nsamples * tsamp:.2f} s")
    print()

    chunks = [
        ("start",    0),
        ("onethird", nsamples_total // 3),
        ("middle",   nsamples_total // 2),
    ]

    all_stats = {}

    for label, start in chunks:
        start = min(start, nsamples_total - args.nsamples)
        print(f"Reading chunk '{label}' at sample {start} "
              f"({start * tsamp:.1f} s)...")
        data = read_chunk(fil_path, start, args.nsamples)

        plot_dynspec(freqs, data, args.flo, args.fhi,
                     start, tsamp, label,
                     outdir / f"dynspec_{label}.png")

        stats = compute_channel_stats(freqs, data, args.flo, args.fhi)
        all_stats[label] = stats

    print("\nWriting statistics:")
    plot_rfi_ratio(all_stats, args.flo, args.fhi,
                   outdir / "dynspec_rfi_ratio.png")
    write_stats_txt(all_stats, args.flo, args.fhi,
                    outdir / "dynspec_stats.txt")

    print("\nDone.")


if __name__ == "__main__":
    main()