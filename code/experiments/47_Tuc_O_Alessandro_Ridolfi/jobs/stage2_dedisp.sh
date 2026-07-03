#!/bin/bash
#SBATCH --job-name=stage2_dedisp
#SBATCH --output=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/slurm_outputs/slurm-dedisp-%A_%a.out
#SBATCH --error=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/slurm_outputs/slurm-dedisp-%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g
#SBATCH --tmp=40g
#SBATCH --time=02:00:00
#SBATCH --array=0-20   # one task per DM trial; must equal len(dm_list.txt) - 1

# ----------------------------------------------------------------------------
# Stage 2: incoherent dedispersion with RFI mitigation using PulsarX
# dedisperse_all_fil inside Apptainer.
#
# Runs one SLURM array task per DM trial. The DM list is read from
# dm_list.txt (one value per line, comments and blank lines ignored).
# Each task produces a dedispersed time series in PRESTO format:
#
#   stage2_dedisp/DM<XX.XX>/
#       <rootname>_DM<XX.XX>.dat
#       <rootname>_DM<XX.XX>.inf
#
# RFI mitigation is applied inline by dedisperse_all_fil using kadaneF and
# zdot algorithms. No separate filtool cleaning step is required because the
# FBFUSE beamformer pipeline pre-normalises the filterbank before writing it.
# The RFI mitigation here catches any residual intermittent RFI.
#
# NOTE: dedisperse_all_fil does not barycentre. The output is topocentric.
# The Viterbi pipeline operates on topocentric time series.
#
# Array size must equal the number of non-comment lines in DM_LIST minus 1.
# ----------------------------------------------------------------------------

set -euo pipefail

# --- paths (edit if needed) -------------------------------------------------
SETUP="/fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh"
EXP_DIR="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/47Tuc_blind_search_v1"
JOBS_DIR="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs"

INPUT_FIL="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/data/47Tuc_22UL_1of2_L.fil"
PULSARX_SIF="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/data/pulsarx_20241219.sif"
BIND_PATH="/fred/oz022/vdimarco"

DM_LIST="${JOBS_DIR}/dm_list.txt"
ROOTNAME="47Tuc"
NTHREADS=4

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
DM_TAG=$(printf "DM%05.2f" "${DM}")

echo "DM:           ${DM}  (tag: ${DM_TAG})"
echo "Input:        ${INPUT_FIL}"

# --- check input exists -----------------------------------------------------
if [ ! -f "${INPUT_FIL}" ]; then
    echo "ERROR: input filterbank not found: ${INPUT_FIL}" >&2
    exit 1
fi

# --- output directory -------------------------------------------------------
OUT_DIR="${EXP_DIR}/stage2_dedisp/${DM_TAG}"
mkdir -p "${OUT_DIR}"

echo "Output dir:   ${OUT_DIR}"

# --- environment ------------------------------------------------------------
set +e
source "${SETUP}"
set -e

module load apptainer

# --- dedisperse -------------------------------------------------------------
echo ""
echo "=== Dedispersion at DM=${DM} with RFI mitigation ==="
cd "${OUT_DIR}"

apptainer exec -B "${BIND_PATH}" "${PULSARX_SIF}" \
    dedisperse_all_fil \
        --dms  "${DM}" \
        --ddm  0.1 \
        --ndm  1 \
        --rfi  kadaneF 8 4 zdot \
        --fillPatch rand \
        --threads "${NTHREADS}" \
        --rootname "${ROOTNAME}" \
        -v \
        --format presto \
        -f "${INPUT_FIL}" \
    2>&1 | tee "${OUT_DIR}/dedisperse.log"

echo ""
echo "Outputs:"
ls -lh "${OUT_DIR}/${ROOTNAME}"*.dat "${OUT_DIR}/${ROOTNAME}"*.inf 2>/dev/null || \
    echo "WARNING: expected output files not found in ${OUT_DIR}"

echo ""
echo "=== Stage 2 task ${IDX} (${DM_TAG}) complete ==="
echo "End (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"