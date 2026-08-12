#!/bin/bash
#SBATCH --job-name=cand153_refold
#SBATCH --output=cand153_refold_%a.out
#SBATCH --error=cand153_refold_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-5

# ============================================================================
# Refold candidate 153 (433.32 Hz) with prepfold acceleration optimisation
# search at a fan of trial DMs, per Emma's suggestion (Slack, 20 Jul).
# Seeds: topocentric P and Pdot from the original DM 0 fold.
# ============================================================================

EXP_DIR="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/47Tuc_blind_search_v1"
SETUP_SCRIPT="/fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh"

CLEANED_FIL="${EXP_DIR}/stage1_clean/47Tuc_22UL_1of2_L_RFIcleaned_01.fil"
OUT_DIR="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/cand_153/cand153_refold"

# Seed parameters (topocentric, from the original prepfold DM 0 fold)
P_SEED="0.00230768641"        # s
PD_SEED="-1.1432e-12"         # s/s

# DM trial fan
DM_TRIALS=(0.0 5.0 10.0 15.0 20.0 25.0)
DM="${DM_TRIALS[${SLURM_ARRAY_TASK_ID}]}"

# ============================================================================
# Provenance header
# ============================================================================
echo "==========================================================================="
echo "Candidate 153 refold, DM trial ${DM}"
echo "Start (UTC):      $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Hostname:         $(hostname)"
echo "Slurm job ID:     ${SLURM_JOB_ID:-N/A}, array task ${SLURM_ARRAY_TASK_ID}"
echo "Input filterbank: ${CLEANED_FIL}"
echo "Seed P (s):       ${P_SEED}"
echo "Seed Pdot (s/s):  ${PD_SEED}"
echo "==========================================================================="

# Environment (guard against Lmod version-swap non-zero exits)
set +e
source "${SETUP_SCRIPT}"
set -e

if [ ! -f "${CLEANED_FIL}" ]; then
    echo "ERROR: cleaned filterbank not found at ${CLEANED_FIL}" >&2
    exit 1
fi

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

# ============================================================================
# prepfold with P, Pdot and DM search enabled (search is on by default)
# ============================================================================
prepfold \
    -p "${P_SEED}" \
    -pd "${PD_SEED}" \
    -dm "${DM}" \
    -npfact 2 \
    -ndmfact 3 \
    -npart 64 \
    -n 64 \
    -fine \
    -noxwin \
    -o "cand153_dm${DM}" \
    "${CLEANED_FIL}"

echo "Done. Best-fit summary:"
grep -A 20 "Best" "cand153_dm${DM}"*.bestprof 2>/dev/null || true
echo "End (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"