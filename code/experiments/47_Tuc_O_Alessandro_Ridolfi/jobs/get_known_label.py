#!/usr/bin/env python
"""
get_known_label.py

Look up the classification label for one candidate (by its dedup-CSV row
index) from classified.csv, and print a short legend label. Prints an empty
line for unclassified candidates (RFI or genuinely new), leaving their
legend unannotated.

    known_main      ->  "<pulsar> (fundamental)"
    known_harmonic  ->  "2x <pulsar> harmonic"      (integer ratio n/1)
                    ->  "3/2 x <pulsar> harmonic"   (fractional ratio)
    unclassified    ->  ""
"""
import argparse, csv

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--idx", type=int, required=True)
    args = ap.parse_args()

    label = ""
    try:
        with open(args.csv) as f:
            for r in csv.DictReader(f):
                if int(r["cand_index"]) == args.idx:
                    cls = r["classification"]
                    par = r["parent_pulsar"]
                    rt = r["harmonic_ratio"]
                    if cls == "known_main":
                        label = f"{par} (fundamental)"
                    elif cls == "known_harmonic":
                        num, den = rt.split("/")
                        if den == "1":
                            label = f"{num}x {par} harmonic"
                        else:
                            label = f"{rt} x {par} harmonic"
                    break
    except FileNotFoundError:
        pass
    print(label)

if __name__ == "__main__":
    main()