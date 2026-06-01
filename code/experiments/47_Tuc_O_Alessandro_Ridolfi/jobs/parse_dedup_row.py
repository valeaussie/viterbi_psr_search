#!/usr/bin/env python
"""
parse_dedup_row.py

Parse one row from candidates_dedup.csv by zero-based row index (excluding
header) and print the fields needed by stage3b_fit.sh as shell variable
assignments, one per line.

The CSV has a structural quirk: the `nt_values` column contains
comma-separated integers (e.g. "16,32,64,128,256") which are written
unquoted, producing a variable number of raw CSV columns. Naive CSV
parsing therefore misidentifies Nt_best, subband_f0, threshold, and
known_match. This script parses correctly by fixing the first three and
last four columns and treating everything in between as nt_values.

Column layout (0-indexed from left / right):
    0:    peak_freq_hz
    1:    peak_loglike
    2:    multiplicity
    3..-5: nt_values (variable, joined back with commas)
    -4:   Nt_best
    -3:   subband_f0
    -2:   threshold
    -1:   known_match

Usage
-----
    python parse_dedup_row.py --csv <candidates_dedup.csv> --idx <N>

Output (stdout, one per line, suitable for eval in bash):
    PEAK_FREQ=<value>
    PEAK_LL=<value>
    MULT=<value>
    NT_VALS=<value>
    NT_BEST=<value>
    SUB_F0=<value>
    KNOWN=<value>

Exit code 1 on any error.
"""

import argparse
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True,
                    help="Path to candidates_dedup.csv")
    ap.add_argument("--idx", type=int, required=True,
                    help="Zero-based row index (excluding header)")
    args = ap.parse_args()

    try:
        with open(args.csv) as fh:
            lines = [l.rstrip("\n") for l in fh]
    except OSError as exc:
        print(f"ERROR: cannot open {args.csv}: {exc}", file=sys.stderr)
        sys.exit(1)

    # Skip header
    data_lines = [l for l in lines[1:] if l.strip()]

    if args.idx >= len(data_lines):
        print(f"ERROR: index {args.idx} >= number of data rows "
              f"{len(data_lines)}", file=sys.stderr)
        sys.exit(1)

    line = data_lines[args.idx]
    parts = line.split(",")

    # Need at least 7 parts: 3 fixed left + at least 1 nt_value + 4 fixed right
    if len(parts) < 7:
        print(f"ERROR: too few columns ({len(parts)}) in row {args.idx}: "
              f"{line!r}", file=sys.stderr)
        sys.exit(1)

    peak_freq  = parts[0].strip()
    peak_ll    = parts[1].strip()
    mult       = parts[2].strip()
    nt_values  = ",".join(p.strip() for p in parts[3:-4])
    nt_best    = parts[-4].strip()
    subband_f0 = parts[-3].strip()
    # parts[-2] is threshold — not needed by the fit script
    known      = parts[-1].strip()

    # Print as shell variable assignments
    print(f"PEAK_FREQ={peak_freq}")
    print(f"PEAK_LL={peak_ll}")
    print(f"MULT={mult}")
    print(f"NT_VALS={nt_values}")
    print(f"NT_BEST={nt_best}")
    print(f"SUB_F0={subband_f0}")
    print(f"KNOWN={known}")


if __name__ == "__main__":
    main()