#!/usr/bin/env python3
"""Normalised Viterbi significance Z for the N_T grid of 47 Tuc O.

For each run, Z = (L_peak - m) / (1.4826 * s), where
  L_peak = max of the per-bin loglikes (delta[:, -1]),
  m      = median of the floor (all bins except the peak window),
  s      = MAD of the floor = median(|L - m|),
1.4826 rescales the MAD to a Gaussian-equivalent sigma.

Reads the *_loglikes.dat files (two columns: freq_Hz, loglike),
prints a table, and writes a CSV and a text table.
"""
import os
import csv
import numpy as np

# ------------------------------------------------------------------ config
DIRPATH   = "."                                  # where the *_loglikes.dat live
PREFIX    = "47TucO_redetect_Nt"                 # file = PREFIX + str(NT) + SUFFIX
SUFFIX    = "_loglikes.dat"
MASK_BINS = 20                                    # half-width of peak mask, in bins
OUT_CSV   = "norm_viterbi_significance.csv"
OUT_TXT   = "norm_viterbi_significance.txt"
RELIABLE_NT_MAX = 128                             # runs with N_T > this: Z unreliable

# grid: N_T and the exact T_coh (= Tsft) used for each run
NTS   = [16, 20, 24, 28, 30, 31, 32, 40, 48, 56, 64, 128, 256, 512, 1024]
TSFT  = [899.713885, 719.771093, 599.809282, 514.122187, 479.847395,
         464.368420, 449.856943, 359.885508, 299.904603, 257.061055,
         224.928433, 112.464217, 56.232108, 28.116016, 14.058008]

# fold S/N you already collated (edit or replace with a live read if needed)
FOLD_SNR = {16: 4.117, 20: 3.276, 24: 4.207, 28: 10.430, 30: 7.139,
            31: 7.580, 32: 7.461, 40: 4.429, 48: 3.685, 56: 4.656,
            64: 7.698, 128: 4.572, 256: 2.987, 512: 3.753, 1024: 3.310}
# ------------------------------------------------------------------

def compute_Z(freq, L, mask_bins=MASK_BINS):
    """Return (Z, peak_freq, L_peak, m, s) from a per-bin loglikes array."""
    ipk = int(np.argmax(L))
    lo  = max(0, ipk - mask_bins)
    hi  = min(len(L), ipk + mask_bins + 1)
    floor = np.concatenate([L[:lo], L[hi:]])
    m = np.median(floor)
    s = np.median(np.abs(floor - m))          # MAD
    if s <= 0:
        return np.nan, freq[ipk], L[ipk], m, s
    Z = (L[ipk] - m) / (1.4826 * s)
    return Z, freq[ipk], L[ipk], m, s

rows = []                                        # collect for file output
hdr = f"{'N_T':>5} {'T_coh(s)':>9} {'peak_f(Hz)':>12} {'L_peak':>10} " \
      f"{'Z(sigma)':>10} {'foldSNR':>8}  note"
lines = [hdr, "-" * len(hdr)]

for NT, tc in zip(NTS, TSFT):
    fold = FOLD_SNR.get(NT, float('nan'))
    path = os.path.join(DIRPATH, f"{PREFIX}{NT}{SUFFIX}")
    if not os.path.isfile(path):
        lines.append(f"{NT:>5} {tc:>9.1f} {'--':>12} {'--':>10} "
                     f"{'MISSING':>10} {fold:>8.2f}  file not found")
        rows.append([NT, tc, "", "", "", fold, False, "file not found"])
        continue
    d = np.loadtxt(path)                          # '#' header skipped automatically
    freq, L = d[:, 0], d[:, 1]
    Z, fpk, Lpk, m, s = compute_Z(freq, L)
    reliable = NT <= RELIABLE_NT_MAX
    note = "" if reliable else "broad mound, Z unreliable"
    lines.append(f"{NT:>5} {tc:>9.1f} {fpk:>12.5f} {Lpk:>10.1f} "
                 f"{Z:>10.2f} {fold:>8.2f}  {note}")
    rows.append([NT, tc, f"{fpk:.6f}", f"{Lpk:.4f}", f"{Z:.4f}",
                 fold, reliable, note])

# ---- print to screen ----
print("\n".join(lines))

# ---- write text table ----
with open(OUT_TXT, "w") as f:
    f.write("\n".join(lines) + "\n")

# ---- write CSV ----
with open(OUT_CSV, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["N_T", "T_coh_s", "peak_freq_Hz", "L_peak",
                "Z_sigma", "fold_snr", "reliable", "note"])
    w.writerows(rows)

print(f"\nWrote {OUT_TXT} and {OUT_CSV}")