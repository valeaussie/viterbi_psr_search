#!/usr/bin/env python
"""
classify_known.py

Classify each blind-search candidate in candidates_dedup.csv against the
known 47 Tuc pulsar catalogue as one of:
    known_main       f_cand ~ a known pulsar F0 (the pulsar itself)
    known_harmonic   f_cand ~ (n/m) * F0 of a known pulsar
    unclassified     no acceptable match (RFI, or genuinely new)

A match is accepted only if BOTH:
    * |f_cand - (n/m)*F0| <= --tol-hz
    * |dm_cand - DM_pulsar| <= --dm-tol            (the DM gate)

Ratios allowed: integers 1..--max-int, plus fractions p/q with
p,q <= --max-frac-pq (drops spurious high-order ratios like 7/8, 7/4).
Winner among accepted matches: lowest harmonic order first, then smallest
frequency residual, then smallest DM residual.

Candidate index matches the dedup-CSV row order (0-based), i.e. the
cand_NNN numbering Stage 3b assigns.
"""
import argparse, csv
from fractions import Fraction

def load_known(path):
    import yaml
    with open(path) as f:
        d = yaml.safe_load(f)
    out = []
    for p in d.get("pulsars", []):
        if p.get("F0_hz") is not None:
            out.append((p["name"], float(p["F0_hz"]),
                        float(p["DM_pccm3"]) if p.get("DM_pccm3") is not None else None))
    return out

def build_ratios(max_int, max_frac_pq):
    rs = set(Fraction(n, 1) for n in range(1, max_int + 1))
    for p in range(1, max_frac_pq + 1):
        for q in range(1, max_frac_pq + 1):
            rs.add(Fraction(p, q))
    return sorted(rs)

def order(fr):
    return fr.numerator + fr.denominator

def classify(f, dm, known, ratios, tol, dm_tol):
    matches = []
    for name, f0, pdm in known:
        if pdm is not None and abs(dm - pdm) > dm_tol:
            continue
        for fr in ratios:
            pred = float(fr) * f0
            dfreq = abs(f - pred)
            if dfreq <= tol:
                ddm = abs(dm - pdm) if pdm is not None else 999.0
                matches.append((order(fr), dfreq, ddm, name, fr, f0, pdm, pred))
    if not matches:
        return None
    matches.sort(key=lambda t: (t[0], t[1], t[2]))
    best = matches[0]
    distinct = len({m[3] for m in matches if m[0] == best[0]})
    return best, distinct, len(matches)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cands", required=True)
    ap.add_argument("--known-yaml", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tol-hz", type=float, default=0.15)
    ap.add_argument("--dm-tol", type=float, default=0.3)
    ap.add_argument("--max-int", type=int, default=8)
    ap.add_argument("--max-frac-pq", type=int, default=3)
    args = ap.parse_args()

    known = load_known(args.known_yaml)
    ratios = build_ratios(args.max_int, args.max_frac_pq)

    with open(args.cands) as f:
        rows = list(csv.DictReader(f))

    out_rows = []
    n_main = n_harm = n_uncl = 0
    for i, r in enumerate(rows):
        f = float(r["peak_freq_hz"]); dm = float(r["dm_best"])
        res = classify(f, dm, known, ratios, args.tol_hz, args.dm_tol)
        if res is None:
            cls, parent, ratio_s, pf0, pdm, pred, dfr, ddm, namb = \
                "unclassified", "", "", "", "", "", "", "", 0
            n_uncl += 1
        else:
            (best, distinct, nmatch) = res
            _, dfreq, ddm, name, fr, f0, pdm_, pred = best
            parent = name; pf0 = f0; pdm = pdm_
            dfr = dfreq; namb = distinct
            if fr == 1:
                cls = "known_main"; ratio_s = "1"; n_main += 1
            else:
                cls = "known_harmonic"; ratio_s = f"{fr.numerator}/{fr.denominator}"; n_harm += 1
        out_rows.append({
            "cand_index": i,
            "peak_freq_hz": f"{f:.6f}",
            "dm_best": f"{dm:.3f}",
            "classification": cls,
            "parent_pulsar": parent,
            "harmonic_ratio": ratio_s,
            "parent_f0_hz": f"{pf0:.6f}" if pf0 != "" else "",
            "parent_dm": f"{pdm:.4f}" if pdm not in ("", None) else "",
            "predicted_freq_hz": f"{pred:.6f}" if pred != "" else "",
            "freq_residual_hz": f"{dfr:.4f}" if dfr != "" else "",
            "dm_residual": f"{ddm:.3f}" if ddm not in ("", None) else "",
            "n_ambiguous": namb,
        })

    cols = ["cand_index","peak_freq_hz","dm_best","classification","parent_pulsar",
            "harmonic_ratio","parent_f0_hz","parent_dm","predicted_freq_hz",
            "freq_residual_hz","dm_residual","n_ambiguous"]
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
        for r in out_rows: w.writerow(r)

    print(f"candidates            : {len(rows)}")
    print(f"  known_main          : {n_main}")
    print(f"  known_harmonic      : {n_harm}")
    print(f"  unclassified (RFI+?) : {n_uncl}")
    print(f"settings: tol={args.tol_hz} Hz, dm_tol={args.dm_tol}, "
          f"int<= {args.max_int}, frac p,q<= {args.max_frac_pq}")
    print(f"\nwrote {args.out}")

    from collections import Counter
    c_main = Counter(r["parent_pulsar"] for r in out_rows if r["classification"]=="known_main")
    c_harm = Counter((r["parent_pulsar"], r["harmonic_ratio"]) for r in out_rows if r["classification"]=="known_harmonic")
    print("\nknown_main by pulsar:")
    for k,v in sorted(c_main.items()): print(f"  {k:14s} {v}")
    print("\nknown_harmonic by (pulsar, ratio):")
    for (k,ratio),v in sorted(c_harm.items()): print(f"  {k:14s} {ratio:>5s}  x{v}")
    amb = [r for r in out_rows if r["n_ambiguous"]>1]
    if amb:
        print(f"\nAMBIGUOUS (more than one pulsar at best order): {len(amb)}")

if __name__ == "__main__":
    main()