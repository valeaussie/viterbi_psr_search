#!/bin/bash
#SBATCH --job-name=stage2_dedisp
#SBATCH --output=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/slurm_outputs/slurm-dedisp-%A_%a.out
#SBATCH --error=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/slurm_outputs/slurm-dedisp-%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16g
#SBATCH --tmp=40g
#SBATCH --time=02:00:00
#SBATCH --array=0-20   # one task per DM trial; must equal len(dm_list.txt) - 1

# ----------------------------------------------------------------------------
# Stage 2: incoherent dedispersion + barycentring with PRESTO prepdata.
#
# Runs one SLURM array task per DM trial.  The DM list is read from
# dm_list.txt (one value per line, comments and blank lines ignored).
# Each task produces a barycentred dedispersed time series:
#
#   stage2_dedisp/DM<XX.XX>_bary/
#       <rootname>_DM<XX.XX>_bary.dat
#       <rootname>_DM<XX.XX>_bary.inf
#
# NOTE: the filterbank channels are already coherently dedispersed within
# each channel at the cluster DM at recording time.  prepdata performs
# the incoherent dedispersion across channels and sums them into a 1-D
# float32 .dat time series.
#
# Array size must equal the number of non-comment lines in DM_LIST minus 1.
# ----------------------------------------------------------------------------

# --- paths (edit if needed) -------------------------------------------------
SETUP="/fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh"
EXP_DIR="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/47Tuc_blind_search_v1"
JOBS_DIR="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs"

CLEANED_FIL="${EXP_DIR}/stage1_clean/47Tuc_22UL_1of2_L_RFIcleaned_01.fil"
DM_LIST="${JOBS_DIR}/dm_list.txt"

# Root name used as prepdata -o prefix (no DM tag; that is appended below)
ROOTNAME="47Tuc"

# ----------------------------------------------------------------------------
IDX=${SLURM_ARRAY_TASK_ID}

echo "Start (UTC):  $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Host:         $(hostname)"
echo "Array job:    ${SLURM_ARRAY_JOB_ID:-N/A} task ${IDX}"

# --- read DM list (strip comments and blank lines) --------------------------
mapfile -t DM_VALS < <(grep -v '^\s*#' "${DM_LIST}" | grep -v '^\s*$')
N_DM=${#DM_VALS[@]}

if [ "${IDX}" -ge "${N_DM}" ]; then
    echo "Array index ${IDX} >= N_DM ${N_DM}; nothing to do."
    exit 0
fi

DM="${DM_VALS[${IDX}]}"

# Format DM tag with two decimal places, e.g. 24.36 -> DM24.36
DM_TAG=$(printf "DM%05.2f" "${DM}")
echo "DM:           ${DM}  (tag: ${DM_TAG})"
echo "Input:        ${CLEANED_FIL}"

# --- environment ------------------------------------------------------------
set +e
source "${SETUP}"
set -e
set -uo pipefail

if ! command -v prepdata >/dev/null 2>&1; then
    echo "ERROR: prepdata not found on PATH after sourcing setup." >&2
    exit 1
fi

if [ ! -f "${CLEANED_FIL}" ]; then
    echo "ERROR: cleaned filterbank not found: ${CLEANED_FIL}" >&2
    exit 1
fi

# --- output directory -------------------------------------------------------
BARY_DIR="${EXP_DIR}/stage2_dedisp/${DM_TAG}_bary"
mkdir -p "${BARY_DIR}"

OUT_ROOT="${ROOTNAME}_${DM_TAG}_bary"

# --- barycentred dedispersion -----------------------------------------------
echo ""
echo "=== Barycentred dedispersion at DM=${DM} ==="
cd "${BARY_DIR}"
prepdata \
    -dm "${DM}" \
    -o "${OUT_ROOT}" \
    "${CLEANED_FIL}" \
    2>&1 | tee "${BARY_DIR}/prepdata.log"

echo ""
echo "Outputs:"
ls -lh "${BARY_DIR}/${OUT_ROOT}.dat" "${BARY_DIR}/${OUT_ROOT}.inf"

echo ""
echo "=== Stage 2 task ${IDX} (${DM_TAG}) complete ==="
echo "End (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Output:    ${BARY_DIR}/${OUT_ROOT}.dat"