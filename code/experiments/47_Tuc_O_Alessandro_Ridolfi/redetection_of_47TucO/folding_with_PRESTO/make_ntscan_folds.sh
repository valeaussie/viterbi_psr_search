#!/bin/bash
# make_ntscan_folds.sh
#
# For every Viterbi run of the coherent-timescale grid (Table 2), i.e. every
# redetection_of_47TucO/47TucO_redetect_Nt<N>_fit_summary.dat, writes the
# polynomial-fit polyco file and prepfold command with make_polyco_fold.py.
# Variants such as Nt28_poly_3f are skipped.
#
# Outputs in folding_with_PRESTO/ntscan/:
#   47TucO_Nt<N>_poly.polycos, 47TucO_Nt<N>_polyco_fold_params.{log,json}
#   prepfold_commands_ntscan.txt   one prepfold command per run, for SLURM

set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
REDETECT=$(dirname "$HERE")
OUT="$HERE/ntscan"
mkdir -p "$OUT"
ALL="$OUT/prepfold_commands_ntscan.txt"
: > "$ALL"

for f in "$REDETECT"/47TucO_redetect_Nt*_fit_summary.dat; do
    b=$(basename "$f")
    if [[ $b =~ ^47TucO_redetect_Nt([0-9]+)_fit_summary\.dat$ ]]; then
        N=${BASH_REMATCH[1]}
    else
        echo "Skipping $b (not a plain Nt run)"
        continue
    fi
    python3 "$HERE/make_polyco_fold.py" --fit "$f" --name "47TucO_Nt${N}" \
        --outdir "$OUT" --folds poly --commands-file "commands_Nt${N}.txt" > /dev/null
    cat "$OUT/commands_Nt${N}.txt" >> "$ALL"
    nfail=$(grep "^CHECK" "$OUT/47TucO_Nt${N}_polyco_fold_params.log" | grep -c "False" || true)
    echo "Nt=${N}: ${nfail} failed checks"
done
echo "$(wc -l < "$ALL") commands written to $ALL"