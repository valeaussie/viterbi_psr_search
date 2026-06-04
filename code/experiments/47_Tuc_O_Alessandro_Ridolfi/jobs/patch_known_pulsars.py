#!/usr/bin/env python3
"""
patch_known_pulsars.py

Patch known_47tuc_pulsars.yaml by adding the 15 new pulsars from
Chen & Risbud (2026, A&A) that are not yet in the ATNF catalogue.

Source for new entries:
    Chen, W., Risbud, D. et al. 2026, A&A, "Fifteen new millisecond
    pulsars in 47 Tucanae", Table 3.

For each new pulsar we record:
    name       JNAME (standardised to J0024-7204xx)
    F0_hz      spin frequency in Hz  (= 1000 / Period_ms from Table 3)
    DM_pccm3   dispersion measure from Table 3
    PB_days    orbital period in days where known, else null
    A1_lts     projected semi-major axis (lt-s) where known, else null
    ECC        eccentricity where known (ai only), else null
    BINARY_type classification string, else null

Fields left as null will not affect the cross-match in aggregate_blind.py,
which only uses name and F0_hz.

Usage
-----
    python patch_known_pulsars.py \\
        --yaml  <path>/config/known_47tuc_pulsars.yaml \\
        --out   <path>/config/known_47tuc_pulsars.yaml   # can overwrite

    Use --dry-run to print what would be added without writing.

After running, inspect the output in VS Code or with:
    python3 -c "import yaml,sys; print(yaml.dump(yaml.safe_load(open('known_47tuc_pulsars.yaml')), default_flow_style=False, sort_keys=False))" | less
"""

import argparse
import copy
import sys

try:
    import yaml
except ImportError:
    sys.exit("ERROR: pyyaml not installed. Run: pip install pyyaml --break-system-packages")

# ---------------------------------------------------------------------------
# New pulsars from Chen & Risbud 2026, Table 3.
# Frequencies computed as 1000 / Period_ms; DMs, Pb, A1 from Table 3.
# ---------------------------------------------------------------------------

NEW_PULSARS = [
    # ae -- He WD binary, Pb = 0.757 d, A1 = 0.741 lt-s
    {
        "name":         "J0024-7204ae",
        "F0_hz":        1000.0 / 3.87,
        "DM_pccm3":     24.34,
        "PB_days":      0.757,
        "A1_lts":       0.741,
        "ECC":          None,
        "BINARY_type":  "He WD",
        "source":       "Chen+2026 Table 3",
    },
    # af -- BW/RB, Pb = 0.0677 d, A1 = 0.0852 lt-s
    {
        "name":         "J0024-7204af",
        "F0_hz":        1000.0 / 2.99,
        "DM_pccm3":     24.34,
        "PB_days":      0.0677,
        "A1_lts":       0.0852,
        "ECC":          None,
        "BINARY_type":  "BW/RB",
        "source":       "Chen+2026 Table 3",
    },
    # ag -- binary, Pb = 1.08 d, A1 = 0.0821 lt-s
    {
        "name":         "J0024-7204ag",
        "F0_hz":        1000.0 / 9.76,
        "DM_pccm3":     24.41,
        "PB_days":      1.08,
        "A1_lts":       0.0821,
        "ECC":          None,
        "BINARY_type":  "binary",
        "source":       "Chen+2026 Table 3",
    },
    # ah -- binary, orbital parameters undetermined
    {
        "name":         "J0024-7204ah",
        "F0_hz":        1000.0 / 3.07,
        "DM_pccm3":     24.36,
        "PB_days":      None,
        "A1_lts":       None,
        "ECC":          None,
        "BINARY_type":  "binary",
        "source":       "Chen+2026 Table 3",
    },
    # ai -- massive companion, eccentric, Pb = 1.65 d, A1 = 5.35 lt-s, e = 0.18
    {
        "name":         "J0024-7204ai",
        "F0_hz":        1000.0 / 13.03,
        "DM_pccm3":     24.47,
        "PB_days":      1.65,
        "A1_lts":       5.35,
        "ECC":          0.18,
        "BINARY_type":  "binary",
        "source":       "Chen+2026 Table 3",
    },
    # aj -- isolated
    {
        "name":         "J0024-7204aj",
        "F0_hz":        1000.0 / 6.36,
        "DM_pccm3":     24.38,
        "PB_days":      None,
        "A1_lts":       None,
        "ECC":          None,
        "BINARY_type":  None,
        "source":       "Chen+2026 Table 3",
    },
    # ak -- binary, orbital parameters undetermined
    {
        "name":         "J0024-7204ak",
        "F0_hz":        1000.0 / 3.52,
        "DM_pccm3":     23.91,
        "PB_days":      None,
        "A1_lts":       None,
        "ECC":          None,
        "BINARY_type":  "binary",
        "source":       "Chen+2026 Table 3",
    },
    # al -- black widow, Pb = 0.157 d, A1 = 0.0206 lt-s
    {
        "name":         "J0024-7204al",
        "F0_hz":        1000.0 / 2.67,
        "DM_pccm3":     24.11,
        "PB_days":      0.157,
        "A1_lts":       0.0206,
        "ECC":          None,
        "BINARY_type":  "BW",
        "source":       "Chen+2026 Table 3",
    },
    # am -- binary, orbital parameters undetermined
    {
        "name":         "J0024-7204am",
        "F0_hz":        1000.0 / 4.16,
        "DM_pccm3":     24.55,
        "PB_days":      None,
        "A1_lts":       None,
        "ECC":          None,
        "BINARY_type":  "binary",
        "source":       "Chen+2026 Table 3",
    },
    # an -- binary, orbital parameters undetermined
    {
        "name":         "J0024-7204an",
        "F0_hz":        1000.0 / 2.61,
        "DM_pccm3":     24.12,
        "PB_days":      None,
        "A1_lts":       None,
        "ECC":          None,
        "BINARY_type":  "binary",
        "source":       "Chen+2026 Table 3",
    },
    # ao -- isolated
    {
        "name":         "J0024-7204ao",
        "F0_hz":        1000.0 / 1.88,
        "DM_pccm3":     23.65,
        "PB_days":      None,
        "A1_lts":       None,
        "ECC":          None,
        "BINARY_type":  None,
        "source":       "Chen+2026 Table 3",
    },
    # ap -- isolated
    {
        "name":         "J0024-7204ap",
        "F0_hz":        1000.0 / 5.11,
        "DM_pccm3":     24.36,
        "PB_days":      None,
        "A1_lts":       None,
        "ECC":          None,
        "BINARY_type":  None,
        "source":       "Chen+2026 Table 3",
    },
    # aq -- binary, orbital parameters undetermined
    {
        "name":         "J0024-7204aq",
        "F0_hz":        1000.0 / 3.04,
        "DM_pccm3":     23.63,
        "PB_days":      None,
        "A1_lts":       None,
        "ECC":          None,
        "BINARY_type":  "binary",
        "source":       "Chen+2026 Table 3",
    },
    # ar -- binary, orbital parameters undetermined
    #       Note: period 9.76 ms same as ag but DM 24.16 vs 24.41 -- two distinct pulsars
    {
        "name":         "J0024-7204ar",
        "F0_hz":        1000.0 / 9.76,
        "DM_pccm3":     24.16,
        "PB_days":      None,
        "A1_lts":       None,
        "ECC":          None,
        "BINARY_type":  "binary",
        "source":       "Chen+2026 Table 3",
    },
    # as -- binary, orbital parameters undetermined
    {
        "name":         "J0024-7204as",
        "F0_hz":        1000.0 / 4.02,
        "DM_pccm3":     24.66,
        "PB_days":      None,
        "A1_lts":       None,
        "ECC":          None,
        "BINARY_type":  "binary",
        "source":       "Chen+2026 Table 3",
    },
]


def format_entry(p):
    """Format a new-pulsar dict into the same schema as ATNF entries."""
    return {
        "name":         p["name"],
        "raj":          None,
        "decj":         None,
        "P0_s":         round(1.0 / p["F0_hz"], 10),
        "F0_hz":        round(p["F0_hz"], 8),
        "DM_pccm3":     p["DM_pccm3"],
        "PB_days":      p["PB_days"],
        "A1_lts":       p["A1_lts"],
        "ECC":          p["ECC"],
        "OM_deg":       None,
        "T0_mjd":       None,
        "BINARY_type":  p["BINARY_type"],
        "Tcoh_max_s":   None,   # not computed here; will be None for entries without A1
        "is_binary":    p["BINARY_type"] is not None,
        "source":       p["source"],
    }


def main():
    ap = argparse.ArgumentParser(
        description="Add 15 new Chen+2026 pulsars to known_47tuc_pulsars.yaml.")
    ap.add_argument("--yaml", required=True,
                    help="Path to existing known_47tuc_pulsars.yaml")
    ap.add_argument("--out", required=True,
                    help="Output path (can be same as --yaml to overwrite)")
    ap.add_argument("--dry-run", action="store_true",
                    help="Print what would be added without writing")
    args = ap.parse_args()

    with open(args.yaml) as fh:
        data = yaml.safe_load(fh)

    existing_names = {p["name"] for p in data.get("pulsars", [])}
    print(f"Existing entries: {len(existing_names)}")

    added = []
    skipped = []
    for p in NEW_PULSARS:
        if p["name"] in existing_names:
            skipped.append(p["name"])
        else:
            entry = format_entry(p)
            added.append(entry)

    print(f"New entries to add : {len(added)}")
    if skipped:
        print(f"Already present (skipped): {skipped}")

    print("\nEntries that will be added:")
    for e in added:
        dm_str = f"{e['DM_pccm3']:.2f}" if e["DM_pccm3"] is not None else "  - "
        pb_str = f"{e['PB_days']:.4f} d" if e["PB_days"] is not None else "     -    "
        print(f"  {e['name']}  F0={e['F0_hz']:>10.5f} Hz  DM={dm_str}  "
              f"Pb={pb_str}  type={e['BINARY_type'] or 'isolated'}")

    if args.dry_run:
        print("\n[dry-run] Nothing written.")
        return

    out_data = copy.deepcopy(data)
    out_data["pulsars"].extend(added)
    out_data["pulsars"].sort(key=lambda p: p["name"])

    # Update metadata if present
    if "metadata" in out_data:
        out_data["metadata"]["n_pulsars"] = len(out_data["pulsars"])
        out_data["metadata"]["notes"] = (
            out_data["metadata"].get("notes", "") +
            " Supplemented with 15 new pulsars from Chen+2026 (A&A)."
        ).strip()

    with open(args.out, "w") as fh:
        yaml.dump(out_data, fh, default_flow_style=False,
                  sort_keys=False, allow_unicode=True)

    print(f"\nWrote {len(out_data['pulsars'])} entries to {args.out}")


if __name__ == "__main__":
    main()