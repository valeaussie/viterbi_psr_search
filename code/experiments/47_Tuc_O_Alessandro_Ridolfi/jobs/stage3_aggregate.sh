#!/bin/bash
#SBATCH --job-name=47tuc-select-v2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
set -eo pipefail
: "${SEARCH_CODE_DIR:?Frozen source directory required}"
: "${SEARCH_OUTPUT_BASE:?New blind output required}"
: "${SEARCH_GRID_MANIFEST:?Frozen grid required}"
: "${SEARCH_RUN_ID:?Frozen run ID required}"
: "${SEARCH_CANDIDATE_DIR:?New candidate output required}"
set +e
source /fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh
setup_status=$?
set -euo pipefail
if (( setup_status != 0 )); then exit "${setup_status}"; fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
JOBS="${SEARCH_CODE_DIR}/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs"
KNOWN="${SEARCH_CODE_DIR}/experiments/47_Tuc_O_Alessandro_Ridolfi/47Tuc_blind_search_v1/config/known_47tuc_pulsars.yaml"
python "${JOBS}/make_search_grid.py" --out "${SEARCH_GRID_MANIFEST}" --validate-results "${SEARCH_OUTPUT_BASE}"
python "${JOBS}/aggregate_blind.py" --blind-dir "${SEARCH_OUTPUT_BASE}" --known-yaml "${KNOWN}" \
    --out-dir "${SEARCH_CANDIDATE_DIR}" --expected-curves "${SEARCH_GRID_SIZE}" \
    --expected-grid "${SEARCH_GRID_MANIFEST}" --run-id "${SEARCH_RUN_ID}" \
    --tail-percentile 95 --false-alarm-rate 0.1 --dedup-tol-hz 0.035
python "${JOBS}/classify_known.py" --cands "${SEARCH_CANDIDATE_DIR}/candidates_timescale_alternatives.csv" \
    --known-yaml "${KNOWN}" --out "${SEARCH_CANDIDATE_DIR}/classified_timescale_alternatives.csv"
python "${JOBS}/flag_rfi.py" --classified-csv "${SEARCH_CANDIDATE_DIR}/classified_timescale_alternatives.csv" \
    --known-yaml "${KNOWN}" --cand-dir "${SEARCH_CANDIDATE_DIR}" --out-csv "${SEARCH_CANDIDATE_DIR}/rfi_flags.csv"
