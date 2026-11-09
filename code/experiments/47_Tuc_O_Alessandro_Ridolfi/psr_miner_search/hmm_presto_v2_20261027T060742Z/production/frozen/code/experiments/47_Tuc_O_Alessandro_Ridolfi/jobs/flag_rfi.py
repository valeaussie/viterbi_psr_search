#!/usr/bin/env python3
"""Frequency heuristics are advisory flags, including at known pulsars. Never exclude or move candidates."""
import argparse
import csv
import os
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


def main():
    ap = argparse.ArgumentParser(description="Flag possible RFI; retain every candidate.")
    ap.add_argument("--classified-csv", required=True)
    ap.add_argument("--known-yaml", required=True)
    ap.add_argument("--cand-dir", required=True)
    ap.add_argument("--out-csv", required=True)
    args = ap.parse_args()
    known = load_known(args.known_yaml)
    harmonics = build_harmonics(known, MAX_HARM, BAND_LO, BAND_HI)
    with open(args.classified_csv, newline="") as fh:
        rows = list(csv.DictReader(fh))
    output=[]
    for r in rows:
        f=float(r["peak_freq_hz"])
        tag,reason=apply_rules(f)
        near=nearest_harmonic(f,harmonics)
        protected=r["classification"] in PROTECTED
        verdict = ("protected_known" if protected else "hold_near_known"
                   if tag and near and abs(near[5]) <= NEAR_KNOWN_TOL else tag or "no_frequency_flag")
        output.append(dict(r, rfi_flag=tag or "", verdict=verdict, reason=reason,
            rfi_policy="flag_only", retained="True", protected_known=protected,
            nearest_known=near[0] if near else "",
            nearest_resid_hz=near[5] if near else "",
            folder_status="present" if os.path.isdir(os.path.join(args.cand_dir,r["candidate_id"])) else "not_yet_created"))
    assert len(output)==len(rows)
    assert [r["candidate_id"] for r in output]==[r["candidate_id"] for r in rows]
    fields=list(dict.fromkeys((list(rows[0]) if rows else ["candidate_id"]) +
        ["rfi_flag","verdict","reason","rfi_policy","retained","protected_known","nearest_known","nearest_resid_hz","folder_status"]))
    with open(args.out_csv,"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(output)
    print(f"Retained {len(output)} candidates; {sum(bool(r['rfi_flag']) for r in output)} frequency flags. No folders moved.")
if __name__ == "__main__": main()
