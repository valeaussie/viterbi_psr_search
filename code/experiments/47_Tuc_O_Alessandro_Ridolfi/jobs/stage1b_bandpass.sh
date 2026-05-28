#!/bin/bash
#SBATCH --job-name=stage1b_bandpass
#SBATCH --output=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/slurm_outputs/slurm-bandpass-%j.out
#SBATCH --error=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/slurm_outputs/slurm-bandpass-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8g
#SBATCH --time=01:00:00

# ----------------------------------------------------------------------------
# Stage 1b: bandpass diagnostic plot + RFI statistics.
# Runs entirely on the host (no container). Sources the project setup script
# so the venv Python (with sigpyproc, numpy, matplotlib, yaml) works.
# ----------------------------------------------------------------------------

# --- paths (edit if needed) -------------------------------------------------
SETUP="/fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh"
POSTPROC="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/stage1_clean_postprocess.py"

INPUT_FIL="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/data/47Tuc_22UL_1of2_L.fil"
STAGE1_DIR="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/47Tuc_blind_search_v1/stage1_clean"
CLEANED_FIL="${STAGE1_DIR}/47Tuc_22UL_1of2_L_RFIcleaned_01.fil"
FILTOOL_LOG="${STAGE1_DIR}/filtool.log"

# --- run --------------------------------------------------------------------
echo "Start (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Host:        $(hostname)"
echo "Job ID:      ${SLURM_JOB_ID:-N/A}"

# Source the environment WITHOUT 'set -e', because Lmod returns a non-zero
# exit code on harmless module warnings (e.g. tempo2 version swap), which
# would otherwise abort the job immediately.
set +e
source "${SETUP}"
set -e

# From here on, fail fast on real errors.
set -uo pipefail

python "${POSTPROC}" \
    --input-fil   "${INPUT_FIL}" \
    --cleaned-fil "${CLEANED_FIL}" \
    --stage1-dir  "${STAGE1_DIR}" \
    --filtool-log "${FILTOOL_LOG}"

echo ""
echo "End (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Outputs:"
echo "  ${STAGE1_DIR}/bandpass.png"
echo "  ${STAGE1_DIR}/rfi_stats.yaml"