#!/usr/bin/env python
"""
flag_rfi.py

Identify RFI among the unclassified candidates of a blind Viterbi search
and optionally move their candidate folders aside.

Only candidates labelled 'unclassified' in classified.csv are ever
considered. Anything labelled known_main or known_harmonic is untouchable.

Three veto rules
----------------
R1  swarm  : SWARM_LO <= f <= SWARM_HI
             One broad interference structure fragmented by the
             deduplication tolerance into many nominal candidates.
R2  comb   : |f - nearest multiple of COMB_SPACING| <= COMB_TOL
             Round-number terrestrial tones. They sit slightly below the
             exact multiple because the search runs on barycentred data.
R3  family : within FAM_TOL of any frequency in FAMILY_CENTRES
             A confirmed terrestrial signal and its fractional harmonics.

Override
--------
A candidate within NEAR_KNOWN_TOL of any n/m harmonic of a catalogue
pulsar is never vetoed. It is held aside for inspection instead. This
stops a real harmonic that happens to fall inside a veto window from
being discarded silently.

Safety check
------------
Before anything is reported or moved, all three rules are applied to the
candidates already labelled known_main or known_harmonic. If any rule
would catch one of them, the script aborts. The known pulsars are the
control group and no rule that damages them is acceptable.

Usage
-----
    python flag_rfi.py \\
        --classified-csv <classified.csv> \\
        --known-yaml     <known_47tuc_pulsars.yaml> \\
        --cand-dir       <candidates/> \\
        --out-csv        <rfi_flags.csv> \\
        [--rfi-dir       <candidates_rfi/>] \\
        [--move]

Without --move the script only reports. Nothing on disk changes.
"""

import argparse
import csv
import os
import shutil
import sys

try:
    import yaml
except ImportError:
    sys.exit("PyYAML is required.")


# ---------------------------------------------------------------------------
# Rule parameters
# ---------------------------------------------------------------------------

SWARM_LO        = 990.0
SWARM_HI        = 1010.0

COMB_SPACING    = 25.0
COMB_TOL        = 0.05

FAMILY_CENTRES  = [216.68, 433.357, 650.03, 866.71]
FAM_TOL         = 0.4

NEAR_KNOWN_TOL  = 0.1
MAX_HARM        = 4

BAND_LO         = 50.0
BAND_HI         = 1059.0

PROTECTED       = ("known_main", "known_harmonic")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_known(yaml_path):
    """Return [(name, F0_hz, DM_pccm3), ...] for catalogue pulsars with an F0."""
    with open(yaml_path) as fh:
        data = yaml.safe_load(fh)
    out = []
    for p in data.get("pulsars", []):
        f0 = p.get("F0_hz")
        if f0 is None:
            continue
        out.append((p.get("name", "?"), float(f0), p.get("DM_pccm3")))
    return out


def build_harmonics(known, max_harm, lo, hi):
    """Return [(name, n, m, freq, dm), ...] for all n/m harmonics inside the band."""
    out = []
    for name, f0, dm in known:
        for n in range(1, max_harm + 1):
            for m in range(1, max_harm + 1):
                f = f0 * n / m
                if lo <= f <= hi:
                    out.append((name, n, m, f, dm))
    return out


def nearest_harmonic(freq, harmonics):
    """Return the (name, n, m, freq, dm, residual) closest to freq, or None."""
    if not harmonics:
        return None
    best = min(harmonics, key=lambda h: abs(freq - h[3]))
    return best + (freq - best[3],)


def apply_rules(freq):
    """Return (tag, reason) if freq is vetoed, else (None, '')."""
    if SWARM_LO <= freq <= SWARM_HI:
        return "rfi_swarm", "inside %.1f to %.1f Hz interference swarm" % (
            SWARM_LO, SWARM_HI)

    k = round(freq / COMB_SPACING) * COMB_SPACING
    if abs(freq - k) <= COMB_TOL:
        return "rfi_comb", "%.4f Hz from %.1f Hz comb line" % (abs(freq - k), k)

    for c in FAMILY_CENTRES:
        if abs(freq - c) <= FAM_TOL:
            return "rfi_family", "%.4f Hz from confirmed family line %.3f Hz" % (
                abs(freq - c), c)

    return None, ""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Flag RFI among unclassified blind-search candidates.")
    ap.add_argument("--classified-csv", required=True)
    ap.add_argument("--known-yaml",     required=True)
    ap.add_argument("--cand-dir",       required=True,
                    help="Directory holding the cand_NNN folders.")
    ap.add_argument("--out-csv",        required=True)
    ap.add_argument("--rfi-dir",        default=None,
                    help="Destination for vetoed folders. "
                         "Default <cand-dir>/../candidates_rfi")
    ap.add_argument("--move", action="store_true",
                    help="Actually move folders. Without this, report only.")
    args = ap.parse_args()

    rfi_dir = args.rfi_dir or os.path.join(
        os.path.dirname(os.path.normpath(args.cand_dir)), "candidates_rfi")

    # ------------------------------------------------------------------
    # Load
    # ------------------------------------------------------------------
    known = load_known(args.known_yaml)
    harmonics = build_harmonics(known, MAX_HARM, BAND_LO, BAND_HI)
    print("Loaded %d catalogue pulsars, %d harmonics in band (n,m <= %d)."
          % (len(known), len(harmonics), MAX_HARM))

    with open(args.classified_csv) as fh:
        rows = list(csv.DictReader(fh))
    print("Read %d candidates from %s" % (len(rows), args.classified_csv))

    for r in rows:
        r["_freq"] = float(r["peak_freq_hz"])
        r["_dm"]   = float(r["dm_best"])
        r["_idx"]  = int(r["cand_index"])

    # ------------------------------------------------------------------
    # Safety check against the control group
    # ------------------------------------------------------------------
    control = [r for r in rows if r["classification"] in PROTECTED]
    violations = []
    for r in control:
        tag, reason = apply_rules(r["_freq"])
        if tag:
            violations.append((r, tag, reason))

    print("\nSafety check: %d protected candidates (%s)."
          % (len(control), ", ".join(PROTECTED)))
    if violations:
        print("\nABORTED. The following rules would catch a known pulsar:")
        for r, tag, reason in violations:
            print("   cand_%03d  %10.4f Hz  %-14s  %s  [%s]"
                  % (r["_idx"], r["_freq"], r["classification"], reason, tag))
        print("\nNo rule that damages the control group is acceptable. "
              "Adjust the parameters and rerun.")
        sys.exit(1)
    print("   PASSED. No rule catches a protected candidate.")

    # ------------------------------------------------------------------
    # Classify the unclassified
    # ------------------------------------------------------------------
    targets = [r for r in rows if r["classification"] == "unclassified"]
    print("\nEvaluating %d unclassified candidates." % len(targets))

    out_rows = []
    counts = {}
    for r in targets:
        f = r["_freq"]
        near = nearest_harmonic(f, harmonics)
        tag, reason = apply_rules(f)

        if tag and near is not None and abs(near[5]) <= NEAR_KNOWN_TOL:
            reason = ("HELD: %.4f Hz from %d/%d x %s, overrides %s"
                      % (abs(near[5]), near[1], near[2], near[0], tag))
            tag = "hold_near_known"
        elif not tag:
            tag = "keep"
            reason = ""

        counts[tag] = counts.get(tag, 0) + 1
        out_rows.append({
            "cand_index":        r["_idx"],
            "peak_freq_hz":      "%.4f" % f,
            "dm_best":           "%.2f" % r["_dm"],
            "verdict":           tag,
            "reason":            reason,
            "nearest_known":     near[0] if near else "",
            "nearest_ratio":     "%d/%d" % (near[1], near[2]) if near else "",
            "nearest_freq_hz":   "%.4f" % near[3] if near else "",
            "nearest_resid_hz":  "%+.4f" % near[5] if near else "",
            "nearest_known_dm":  "%.3f" % near[4] if near and near[4] is not None else "",
        })

    out_rows.sort(key=lambda d: float(d["peak_freq_hz"]))

    # ------------------------------------------------------------------
    # Report
    # ------------------------------------------------------------------
    print("\nVerdict summary")
    for k in sorted(counts):
        print("   %-16s %4d" % (k, counts[k]))

    survivors = [d for d in out_rows if d["verdict"] in ("keep", "hold_near_known")]
    if survivors:
        print("\nSurvivors (%d), sorted by frequency:" % len(survivors))
        print("   %-10s %12s %8s  %-16s %-14s %8s %10s"
              % ("cand", "freq_hz", "dm", "verdict",
                 "nearest_known", "ratio", "resid_hz"))
        for d in survivors:
            print("   cand_%03d   %12s %8s  %-16s %-14s %8s %10s"
                  % (int(d["cand_index"]), d["peak_freq_hz"], d["dm_best"],
                     d["verdict"], d["nearest_known"], d["nearest_ratio"],
                     d["nearest_resid_hz"]))

    cols = ["cand_index", "peak_freq_hz", "dm_best", "verdict", "reason",
            "nearest_known", "nearest_ratio", "nearest_freq_hz",
            "nearest_resid_hz", "nearest_known_dm"]
    with open(args.out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for d in out_rows:
            w.writerow(d)
    print("\nWrote %d rows to %s" % (len(out_rows), args.out_csv))

    # ------------------------------------------------------------------
    # Move
    # ------------------------------------------------------------------
    to_move = [d for d in out_rows if d["verdict"].startswith("rfi_")]

    if not args.move:
        print("\nDRY RUN. %d folders would move to %s" % (len(to_move), rfi_dir))
        print("Rerun with --move to apply.")
        return

    if not os.path.isdir(rfi_dir):
        os.makedirs(rfi_dir)
        print("\nCreated %s" % rfi_dir)

    moved = skipped = 0
    for d in to_move:
        name = "cand_%03d" % int(d["cand_index"])
        src = os.path.join(args.cand_dir, name)
        dst = os.path.join(rfi_dir, name)
        if not os.path.isdir(src):
            print("   MISSING  %s" % src)
            skipped += 1
            continue
        if os.path.exists(dst):
            print("   EXISTS   %s, not overwritten" % dst)
            skipped += 1
            continue
        shutil.move(src, dst)
        moved += 1

    print("\nMoved %d folders to %s" % (moved, rfi_dir))
    if skipped:
        print("Skipped %d. See messages above." % skipped)

    manifest = os.path.join(rfi_dir, "rfi_manifest.csv")
    with open(manifest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for d in to_move:
            w.writerow(d)
    print("Wrote manifest to %s" % manifest)


if __name__ == "__main__":
    main()