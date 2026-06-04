#!/bin/bash
#SBATCH --job-name=stage4_fold
#SBATCH --output=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/slurm_outputs/slurm-fold-%A_%a.out
#SBATCH --error=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/slurm_outputs/slurm-fold-%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4g
#SBATCH --time=00:20:00
#SBATCH --array=0-104

# ----------------------------------------------------------------------------
# Stage 4 (fold): fold each Viterbi candidate with PulsarX (psrfold_fil).
#
# Each candidate is folded twice (polynomial and Kepler fits). Each fold
# runs in its own subdirectory so PulsarX auto-named PNG files do not
# overwrite each other:
#
#   FOLD_DIR/
#     cand_NNN/
#       poly/
#         J0000-00_...ar
#         J0000-00_...png    <-- diagnostic plot
#         J0000-00_...cands
#       kepler/
#         J0000-00_...ar
#         J0000-00_...png
#         J0000-00_...cands
#
# No -o flag is passed to psrfold_fil so PulsarX saves PNGs automatically.
# ----------------------------------------------------------------------------

# --- paths ------------------------------------------------------------------
SETUP="/fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh"
CODE_DIR="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code"
EXP_DIR="${CODE_DIR}/experiments/47_Tuc_O_Alessandro_Ridolfi/47Tuc_blind_search_v1"

PULSARX_SIF="${CODE_DIR}/experiments/47_Tuc_O_Alessandro_Ridolfi/data/pulsarx_20241219.sif"
BIND_PATH="/fred/oz022/vdimarco"

FIL="${CODE_DIR}/experiments/47_Tuc_O_Alessandro_Ridolfi/data/47Tuc_22UL_1of2_L.fil"
CAND_DIR="${EXP_DIR}/stage3_viterbi/candidates"
FOLD_DIR="${EXP_DIR}/stage4_fold"
TEMPLATE="/home/pulsarx/software/PulsarX/include/template/meerkat_fold.template"

NTHREADS=2
CLFD=2.0
NBINS=256
NSUBINT=64
SUBLENGTH=30

# --- map array index to candidate directory ---------------------------------
mapfile -t CAND_DIRS < <(find "${CAND_DIR}" -maxdepth 1 -mindepth 1 -type d | sort)

N_CANDS=${#CAND_DIRS[@]}
IDX=${SLURM_ARRAY_TASK_ID}

if [ "${IDX}" -ge "${N_CANDS}" ]; then
    echo "Array index ${IDX} >= N_CANDS ${N_CANDS}; nothing to do."
    exit 0
fi

THIS_CAND_DIR="${CAND_DIRS[${IDX}]}"
CAND_ID=$(basename "${THIS_CAND_DIR}")

echo "Start (UTC):    $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Host:           $(hostname)"
echo "Array job:      ${SLURM_ARRAY_JOB_ID:-N/A} task ${IDX}"
echo "Candidate:      ${CAND_ID}"

# --- environment ------------------------------------------------------------
set +e
source "${SETUP}"
set -e
set -uo pipefail
ml singularity

# --- locate candfile --------------------------------------------------------
CANDFILE=$(find "${THIS_CAND_DIR}" -maxdepth 1 -name "*_psrfold.candfile" | head -1)
INFOFILE="${CANDFILE}.info"

if [ -z "${CANDFILE}" ]; then
    echo "ERROR: no *_psrfold.candfile found in ${THIS_CAND_DIR}"
    exit 1
fi

if [ ! -f "${INFOFILE}" ]; then
    echo "ERROR: sidecar info file not found: ${INFOFILE}"
    exit 1
fi

echo "Candfile:  ${CANDFILE}"

# --- parse pepoch and F2 from sidecar ---------------------------------------
mapfile -t PEPOCHS < <(grep -o -- '--pepoch [0-9.]*' "${INFOFILE}" | awk '{print $2}')
mapfile -t F2S     < <(grep -o -- '--f2 [0-9eE.+-]*' "${INFOFILE}" | awk '{print $2}')

PEPOCH_POLY="${PEPOCHS[0]:-}"
PEPOCH_KEP="${PEPOCHS[1]:-}"
F2_POLY="${F2S[0]:-0.0}"
F2_KEP="${F2S[1]:-0.0}"

if [ -z "${PEPOCH_POLY}" ] || [ -z "${PEPOCH_KEP}" ]; then
    echo "ERROR: could not parse pepoch values from ${INFOFILE}"
    exit 1
fi

echo "pepoch (poly):   ${PEPOCH_POLY}"
echo "pepoch (kepler): ${PEPOCH_KEP}"

# --- parse candfile rows ----------------------------------------------------
mapfile -t CAND_ROWS < <(grep -v '^#' "${CANDFILE}" | grep -v '^[[:space:]]*$')

if [ "${#CAND_ROWS[@]}" -lt 1 ]; then
    echo "ERROR: no data rows found in ${CANDFILE}"
    exit 1
fi

# --- fold function ----------------------------------------------------------
# Each fold runs in its own subdirectory so auto-named PNGs don't overwrite.
fold_candidate() {
    local ROW_ID="$1"
    local ROW="$2"
    local PEPOCH="$3"
    local F2="$4"
    local LABEL="$5"

    local DM F0 F1
    DM=$(echo "${ROW}" | awk '{print $2}')
    F0=$(echo "${ROW}" | awk '{print $4}')
    F1=$(echo "${ROW}" | awk '{print $5}')

    # Each fold gets its own subdirectory
    local SUB_DIR="${FOLD_DIR}/${CAND_ID}/${LABEL}"
    mkdir -p "${SUB_DIR}"
    cd "${SUB_DIR}"

    echo ""
    echo "--- Folding row ${ROW_ID} (${LABEL}) ---"
    echo "  DM=${DM}  F0=${F0}  F1=${F1}  F2=${F2}  pepoch=${PEPOCH}"
    echo "  Output dir: ${SUB_DIR}"

    singularity exec -B "${BIND_PATH}" "${PULSARX_SIF}" psrfold_fil \
        -v \
        -t "${NTHREADS}" \
        --f0     "${F0}" \
        --f1     "${F1}" \
        --f2     "${F2}" \
        --dm     "${DM}" \
        --pepoch "${PEPOCH}" \
        --clfd   "${CLFD}" \
        --rfi kadaneF 8 4 zdot \
        --fillPatch rand \
        -n "${NBINS}" \
        -b "${NSUBINT}" \
        -L "${SUBLENGTH}" \
        --template "${TEMPLATE}" \
        -f "${FIL}"

    # Remove any core dumps
    rm -f "${SUB_DIR}"/core.*

    echo "  Done row ${ROW_ID}."
}

# --- fold poly --------------------------------------------------------------
fold_candidate 0 "${CAND_ROWS[0]}" "${PEPOCH_POLY}" "${F2_POLY}" "poly"

# --- fold kepler (if present) -----------------------------------------------
if [ "${#CAND_ROWS[@]}" -ge 2 ]; then
    fold_candidate 1 "${CAND_ROWS[1]}" "${PEPOCH_KEP}" "${F2_KEP}" "kepler"
fi

echo ""
echo "=== task ${IDX} (${CAND_ID}) complete ==="
echo "End (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"