#!/usr/bin/env python
"""
rank_candidates.py

Filter and rank blind search candidates from candidates_dedup.csv.

Steps
-----
1. Remove candidates matched to a known pulsar (known_match non-empty).
2. Remove candidates that are harmonics of known pulsars (n/m * F0_known
   within --harmonic-tol-hz for integer ratios n/m up to --max-harmonic).
3. Remove candidates that are harmonics of other surviving candidates
   (promotes the one with higher loglike as the fundamental).
4. Score each surviving candidate and rank by score.
5. Write ranked_candidates.csv.

Scoring
-------
Each candidate is scored on four criteria, each normalised to [0,1]:

    dm_proximity   : 1 - |dm_best - cluster_dm| / dm_range_half
                     (1 = exactly at cluster DM, 0 = at edge of DM grid)
    dm_consistency : dm_count / max(dm_count)
                     (how many DM trials recovered it)
    multiplicity   : log(multiplicity) / log(max(multiplicity))
                     (how many subband/Nt combinations recovered it, log-scaled)
    nt_coverage    : len(nt_values.split(';')) / max_possible_nt
                     (fraction of Nt values where it was seen)
    snr_excess     : (peak_loglike - threshold) / max(peak_loglike - threshold)
                     (how far above threshold, normalised)

Final score = weighted sum with weights --w-dm-prox, --w-dm-cons, --w-mult,
--w-nt, --w-snr (defaults 2, 2, 1, 1, 1).

The dm_proximity weight is doubled by default because DM consistency with
the cluster is the single strongest discriminator against RFI (which appears
at DM=0 or random DM) and against noise fluctuations.

Usage
-----
    python rank_candidates.py \\
        --dedup-csv   <candidates_dedup.csv> \\
        --known-yaml  <known_47tuc_pulsars.yaml> \\
        --out-csv     <ranked_candidates.csv> \\
        [--cluster-dm      24.36] \\
        [--dm-range-half   1.2] \\
        [--harmonic-tol-hz 0.1] \\
        [--max-harmonic    6] \\
        [--w-dm-prox  2] \\
        [--w-dm-cons  2] \\
        [--w-mult     1] \\
        [--w-nt       1] \\
        [--w-snr      1]
"""

import argparse
import csv
import math
import os

try:
    import yaml
    HAVE_YAML = True
except ImportError:
    HAVE_YAML = False


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_known_f0s(yaml_path):
    """Return list of (name, F0_hz) for all known pulsars with an F0."""
    if not (HAVE_YAML and yaml_path and os.path.exists(yaml_path)):
        return []
    with open(yaml_path) as fh:
        data = yaml.safe_load(fh)
    out = []
    for p in data.get("pulsars", []):
        f0 = p.get("F0_hz")
        name = p.get("name", "?")
        if f0 is not None:
            out.append((name, float(f0)))
    return out


def is_harmonic_of(freq, known_f0s, tol_hz, max_harmonic):
    """
    Return (True, reason_string) if freq is within tol_hz of any ratio
    n/m * F0 for integers 1 <= n,m <= max_harmonic with n != m or m==1.
    Only ratios that are simple fractions are tested (n and m both <= max_harmonic).
    Returns (False, '') otherwise.
    """
    for name, f0 in known_f0s:
        for n in range(1, max_harmonic + 1):
            for m in range(1, max_harmonic + 1):
                if n == m:
                    continue  # same frequency, already caught by known_match
                harmonic_freq = f0 * n / m
                if abs(freq - harmonic_freq) <= tol_hz:
                    return True, f"{n}/{m} x {name} ({f0:.4f} Hz)"
    return False, ""


def is_harmonic_of_candidate(freq, candidates, tol_hz, max_harmonic):
    """
    Return (True, reason_string) if freq is within tol_hz of any ratio
    n/m * f0 for integers 1 <= n,m <= max_harmonic with n != m,
    for any surviving candidate frequency f0.
    The candidate with higher loglike is taken as the fundamental.
    """
    for c in candidates:
        f0 = c["peak_freq_hz"]
        if abs(freq - f0) < tol_hz:
            continue  # same candidate
        for n in range(1, max_harmonic + 1):
            for m in range(1, max_harmonic + 1):
                if n == m:
                    continue
                harmonic_freq = f0 * n / m
                if abs(freq - harmonic_freq) <= tol_hz:
                    return True, f"{n}/{m} x candidate at {f0:.4f} Hz"
    return False, ""


def nt_count(nt_values_str):
    """Count distinct Nt values from a semicolon-separated string."""
    if not nt_values_str:
        return 0
    return len([x for x in nt_values_str.split(";") if x.strip()])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Filter and rank blind search candidates.")
    ap.add_argument("--dedup-csv",  required=True)
    ap.add_argument("--known-yaml", default=None)
    ap.add_argument("--out-csv",    required=True)

    ap.add_argument("--cluster-dm",      type=float, default=24.36,
                    help="Expected DM of cluster pulsars (default 24.36).")
    ap.add_argument("--dm-range-half",   type=float, default=1.2,
                    help="Half-width of DM search grid (default 1.2 pc/cm^3).")
    ap.add_argument("--harmonic-tol-hz", type=float, default=0.1,
                    help="Frequency tolerance for harmonic check (default 0.1 Hz).")
    ap.add_argument("--max-harmonic",    type=int,   default=6,
                    help="Maximum harmonic ratio n or m to check (default 6).")

    ap.add_argument("--w-dm-prox", type=float, default=2.0)
    ap.add_argument("--w-dm-cons", type=float, default=2.0)
    ap.add_argument("--w-mult",    type=float, default=1.0)
    ap.add_argument("--w-nt",      type=float, default=1.0)
    ap.add_argument("--w-snr",     type=float, default=1.0)
    args = ap.parse_args()

    # ------------------------------------------------------------------
    # Load inputs
    # ------------------------------------------------------------------
    known_f0s = load_known_f0s(args.known_yaml)
    print(f"Loaded {len(known_f0s)} known pulsars.")

    with open(args.dedup_csv) as fh:
        reader = csv.DictReader(fh)
        all_rows = list(reader)
    print(f"Read {len(all_rows)} candidates from {args.dedup_csv}")

    # Cast numeric fields
    for r in all_rows:
        r["peak_freq_hz"]  = float(r["peak_freq_hz"])
        r["peak_loglike"]  = float(r["peak_loglike"])
        r["dm_best"]       = float(r["dm_best"])
        r["dm_count"]      = int(r["dm_count"])
        r["multiplicity"]  = int(r["multiplicity"])
        r["threshold"]     = float(r["threshold"])

    # ------------------------------------------------------------------
    # Step 1: remove known pulsar matches
    # ------------------------------------------------------------------
    after_known = [r for r in all_rows if not r.get("known_match", "").strip()]
    n_known = len(all_rows) - len(after_known)
    print(f"\nStep 1: removed {n_known} known-pulsar matches -> {len(after_known)} remain")

    # ------------------------------------------------------------------
    # Step 2: remove harmonics of known pulsars
    # ------------------------------------------------------------------
    after_known_harm = []
    known_harm_removed = []
    for r in after_known:
        hit, reason = is_harmonic_of(
            r["peak_freq_hz"], known_f0s,
            args.harmonic_tol_hz, args.max_harmonic
        )
        if hit:
            r["reject_reason"] = f"harmonic of known pulsar: {reason}"
            known_harm_removed.append(r)
        else:
            after_known_harm.append(r)

    print(f"Step 2: removed {len(known_harm_removed)} harmonics of known pulsars "
          f"-> {len(after_known_harm)} remain")
    for r in known_harm_removed:
        print(f"  {r['peak_freq_hz']:>10.4f} Hz  {r['reject_reason']}")

    # ------------------------------------------------------------------
    # Step 3: remove inter-candidate harmonics
    # Sort by loglike descending first so the higher-loglike candidate
    # is always kept as the fundamental.
    # ------------------------------------------------------------------
    survivors = sorted(after_known_harm, key=lambda r: r["peak_loglike"], reverse=True)
    after_cand_harm = []
    cand_harm_removed = []
    for r in survivors:
        hit, reason = is_harmonic_of_candidate(
            r["peak_freq_hz"], after_cand_harm,
            args.harmonic_tol_hz, args.max_harmonic
        )
        if hit:
            r["reject_reason"] = f"harmonic of surviving candidate: {reason}"
            cand_harm_removed.append(r)
        else:
            after_cand_harm.append(r)

    print(f"Step 3: removed {len(cand_harm_removed)} inter-candidate harmonics "
          f"-> {len(after_cand_harm)} remain")
    for r in cand_harm_removed:
        print(f"  {r['peak_freq_hz']:>10.4f} Hz  {r['reject_reason']}")

    # ------------------------------------------------------------------
    # Step 4: score and rank survivors
    # ------------------------------------------------------------------
    if not after_cand_harm:
        print("\nNo candidates survived filtering.")
        return

    max_dm_count    = max(r["dm_count"]    for r in after_cand_harm)
    max_mult        = max(r["multiplicity"] for r in after_cand_harm)
    max_nt          = max(nt_count(r["nt_values"]) for r in after_cand_harm)
    max_snr_excess  = max(r["peak_loglike"] - r["threshold"] for r in after_cand_harm)

    W = (args.w_dm_prox + args.w_dm_cons + args.w_mult + args.w_nt + args.w_snr)

    for r in after_cand_harm:
        # DM proximity to cluster
        dm_offset = abs(r["dm_best"] - args.cluster_dm)
        dm_prox = max(0.0, 1.0 - dm_offset / args.dm_range_half)

        # DM consistency (how many trials recovered it)
        dm_cons = r["dm_count"] / max_dm_count

        # Multiplicity (log-scaled to avoid J dominating everything)
        mult = math.log1p(r["multiplicity"]) / math.log1p(max_mult)

        # Nt coverage
        nt = nt_count(r["nt_values"]) / max_nt if max_nt > 0 else 0.0

        # SNR excess above threshold
        snr = (r["peak_loglike"] - r["threshold"]) / max_snr_excess \
              if max_snr_excess > 0 else 0.0

        score = (args.w_dm_prox * dm_prox
                 + args.w_dm_cons * dm_cons
                 + args.w_mult   * mult
                 + args.w_nt     * nt
                 + args.w_snr    * snr) / W

        r["score"]          = round(score, 4)
        r["s_dm_proximity"] = round(dm_prox, 3)
        r["s_dm_consist"]   = round(dm_cons, 3)
        r["s_multiplicity"] = round(mult,    3)
        r["s_nt_coverage"]  = round(nt,      3)
        r["s_snr_excess"]   = round(snr,     3)

    ranked = sorted(after_cand_harm, key=lambda r: r["score"], reverse=True)

    # ------------------------------------------------------------------
    # Write output CSV
    # ------------------------------------------------------------------
    out_cols = [
        "rank", "peak_freq_hz", "peak_loglike", "threshold",
        "dm_best", "dm_count", "dm_values",
        "multiplicity", "nt_values", "Nt_best", "subband_f0",
        "score",
        "s_dm_proximity", "s_dm_consist", "s_multiplicity",
        "s_nt_coverage", "s_snr_excess",
    ]
    with open(args.out_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=out_cols, extrasaction="ignore")
        writer.writeheader()
        for i, r in enumerate(ranked, 1):
            r["rank"] = i
            writer.writerow(r)

    print(f"\nWrote {len(ranked)} ranked candidates to {args.out_csv}")
    print(f"\nTop 15 candidates:")
    hdr = (f"  {'rank':>4}  {'freq_hz':>12}  {'score':>6}  {'dm_best':>7}  "
           f"{'dm_cnt':>6}  {'mult':>5}  {'nt_values':>15}  {'loglike':>9}")
    print(hdr)
    for r in ranked[:15]:
        print(f"  {r['rank']:>4}  {r['peak_freq_hz']:>12.4f}  {r['score']:>6.3f}  "
              f"{r['dm_best']:>7.2f}  {r['dm_count']:>6}  {r['multiplicity']:>5}  "
              f"{r['nt_values']:>15}  {r['peak_loglike']:>9.1f}")


if __name__ == "__main__":
    main()