#!/bin/bash
#SBATCH --job-name=stage5_collate
#SBATCH --output=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/slurm_outputs/slurm-collate-%j.out
#SBATCH --error=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/slurm_outputs/slurm-collate-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=00:10:00

# ----------------------------------------------------------------------------
# Stage 5 (collate): extract S/N from all folded .ar files and produce a
# ranked candidate CSV.
#
# Sources setup_viterbi_psr.sh FIRST so that psrstat is on PATH before
# Python is invoked. collate_fold_snr.py then calls psrstat directly
# in a single subprocess (all files at once).
# ----------------------------------------------------------------------------

SETUP="/fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh"
CODE_DIR="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code"
EXP_DIR="${CODE_DIR}/experiments/47_Tuc_O_Alessandro_Ridolfi/47Tuc_blind_search_v1"

COLLATE_SCRIPT="${CODE_DIR}/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/collate_fold_snr.py"
FOLD_DIR="${EXP_DIR}/stage4_fold"
CAND_DIR="${EXP_DIR}/stage3_viterbi/candidates"
DEDUP_CSV="${EXP_DIR}/stage3_viterbi/blind_v1/candidates_dedup.csv"
OUT_CSV="${FOLD_DIR}/fold_snr_ranked.csv"

echo "Start (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Host:        $(hostname)"
echo "Job ID:      ${SLURM_JOB_ID:-N/A}"

# Source environment BEFORE calling Python so psrstat is on PATH
set +e
source "${SETUP}"
set -e
set -uo pipefail

echo "psrstat location: $(which psrstat)"

python -u "${COLLATE_SCRIPT}" \
    --fold-dir  "${FOLD_DIR}" \
    --cand-dir  "${CAND_DIR}" \
    --dedup-csv "${DEDUP_CSV}" \
    --setup     "${SETUP}" \
    --out-csv   "${OUT_CSV}"

echo ""
echo "=== collate complete ==="
echo "End (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Ranked CSV: ${OUT_CSV}"