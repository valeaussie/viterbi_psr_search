#!/usr/bin/env python
"""
parse_dedup_row.py

Parse one row from candidates_dedup.csv by zero-based row index (excluding
the header) and print the fields needed by stage3b_fit.sh as shell variable
assignments, one per line, suitable for eval in bash.

Column layout of candidates_dedup.csv (written by aggregate_blind.py):

    0:  peak_freq_hz
    1:  peak_loglike
    2:  dm_best
    3:  dm_count
    4:  dm_values       (semicolon-separated, single CSV field)
    5:  multiplicity
    6:  nt_values       (semicolon-separated, single CSV field)
    7:  Nt_best
    8:  subband_f0
    9:  threshold
    10: known_match

dm_values and nt_values use semicolons as internal separators so the CSV
remains exactly 11 columns per row with no ambiguity.

Usage
-----
    python parse_dedup_row.py --csv candidates_dedup.csv --idx <N>

Output (stdout, one per line):
    PEAK_FREQ=<value>
    PEAK_LL=<value>
    DM_BEST=<value>
    DM_COUNT=<value>
    DM_VALUES=<value>
    MULT=<value>
    NT_VALS=<value>
    NT_BEST=<value>
    SUB_F0=<value>
    KNOWN=<value>

Exit code 1 on any error.
"""

import argparse
import sys


def parse_row(line):
    """
    Parse one data line from candidates_dedup.csv.

    Returns a dict with keys:
        peak_freq_hz, peak_loglike, dm_best, dm_count, dm_values,
        multiplicity, nt_values, Nt_best, subband_f0, threshold, known_match

    Raises ValueError with a description if parsing fails.
    """
    parts = line.split(",")

    if len(parts) != 11:
        raise ValueError(
            f"expected 11 columns, got {len(parts)}. "
            f"Check that dm_values and nt_values use semicolons as separators.")

    return {
        "peak_freq_hz": parts[0].strip(),
        "peak_loglike": parts[1].strip(),
        "dm_best":      parts[2].strip(),
        "dm_count":     parts[3].strip(),
        "dm_values":    parts[4].strip(),
        "multiplicity": parts[5].strip(),
        "nt_values":    parts[6].strip(),
        "Nt_best":      parts[7].strip(),
        "subband_f0":   parts[8].strip(),
        "threshold":    parts[9].strip(),
        "known_match":  parts[10].strip(),
    }


def main():
    ap = argparse.ArgumentParser(
        description="Parse one row from candidates_dedup.csv for use in bash.")
    ap.add_argument("--csv", required=True,
                    help="Path to candidates_dedup.csv")
    ap.add_argument("--idx", type=int, required=True,
                    help="Zero-based row index (excluding header)")
    args = ap.parse_args()

    try:
        with open(args.csv) as fh:
            lines = [line.rstrip("\n") for line in fh]
    except OSError as exc:
        print(f"ERROR: cannot open {args.csv}: {exc}", file=sys.stderr)
        sys.exit(1)

    data_lines = [line for line in lines[1:] if line.strip()]

    if args.idx >= len(data_lines):
        print(f"ERROR: index {args.idx} >= number of data rows "
              f"({len(data_lines)})", file=sys.stderr)
        sys.exit(1)

    line = data_lines[args.idx]

    try:
        fields = parse_row(line)
    except ValueError as exc:
        print(f"ERROR: failed to parse row {args.idx}: {exc}", file=sys.stderr)
        sys.exit(1)

    print(f"PEAK_FREQ={fields['peak_freq_hz']}")
    print(f"PEAK_LL={fields['peak_loglike']}")
    print(f"DM_BEST={fields['dm_best']}")
    print(f"DM_COUNT={fields['dm_count']}")
    print(f"DM_VALUES={fields['dm_values']}")
    print(f"MULT={fields['multiplicity']}")
    print(f"NT_VALS={fields['nt_values']}")
    print(f"NT_BEST={fields['Nt_best']}")
    print(f"SUB_F0={fields['subband_f0']}")
    print(f"KNOWN={fields['known_match']}")


if __name__ == "__main__":
    main()