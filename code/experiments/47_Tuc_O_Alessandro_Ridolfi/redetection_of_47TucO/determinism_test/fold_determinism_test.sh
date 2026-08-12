#!/bin/bash
#SBATCH --job-name=det_test
#SBATCH --output=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/slurm_outputs/slurm-det-%A_%a.out
#SBATCH --error=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/slurm_outputs/slurm-det-%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4g
#SBATCH --time=00:40:00
#SBATCH --array=0-23

# ----------------------------------------------------------------------------
# Fold determinism test.
#
# Runs the SAME psrfold_fil command many times to test whether the fold is
# reproducible. Two arms, 12 runs each:
#   tasks  0-11 : --fill rand --fillPatch rand   (current default behaviour)
#   tasks 12-23 : --fill mean --fillPatch mean   (no RNG called)
#
# XLibs seeds its RNG with std::random_device on every run (patch.cpp,
# rfi.cpp), so the 'rand' arm is expected to vary and the 'mean' arm is
# expected to be bit-identical.
#
# Compare afterwards with the checksums (see bottom of this file).
# ----------------------------------------------------------------------------

# --- paths ------------------------------------------------------------------
SETUP="/fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh"
CODE_DIR="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code"
EXP_ROOT="${CODE_DIR}/experiments/47_Tuc_O_Alessandro_Ridolfi"

PULSARX_SIF="${EXP_ROOT}/data/pulsarx_20241219.sif"
BIND_PATH="/fred/oz022/vdimarco"
FIL="${EXP_ROOT}/data/47Tuc_22UL_1of2_L.fil"
TEMPLATE="/home/pulsarx/software/PulsarX/include/template/meerkat_fold.template"

OUT_ROOT="${EXP_ROOT}/redetection_of_47TucO/determinism_test"

# --- fold parameters (f0-only polynomial value, topocentric) ----------------
F0="378.326969588272"
DM="24.356441"
PEPOCH="59391.2119067078"

NTHREADS=2
CLFD=2.0
NBINS=256
NSUBINT=64
SUBLENGTH=30

# --- map array index to arm and repeat number -------------------------------
IDX=${SLURM_ARRAY_TASK_ID}

if [ "${IDX}" -lt 12 ]; then
    ARM="rand"
    REP=${IDX}
    FILL_FLAGS="--fill rand --fillPatch rand"
else
    ARM="mean"
    REP=$(( IDX - 12 ))
    FILL_FLAGS="--fill mean --fillPatch mean"
fi

RUN_DIR="${OUT_ROOT}/${ARM}/rep_$(printf '%02d' ${REP})"

echo "Start (UTC):  $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Host:         $(hostname)"
echo "Array job:    ${SLURM_ARRAY_JOB_ID:-N/A} task ${IDX}"
echo "Arm:          ${ARM}  (repeat ${REP})"
echo "Fill flags:   ${FILL_FLAGS}"
echo "Run dir:      ${RUN_DIR}"

# --- environment ------------------------------------------------------------
set +e
source "${SETUP}"
set -e
set -uo pipefail
ml singularity

mkdir -p "${RUN_DIR}"
cd "${RUN_DIR}"

singularity exec -B "${BIND_PATH}" "${PULSARX_SIF}" psrfold_fil \
    -v \
    -t "${NTHREADS}" \
    --f0     "${F0}" \
    --dm     "${DM}" \
    --pepoch "${PEPOCH}" \
    --clfd   "${CLFD}" \
    --rfi kadaneF 8 4 zdot \
    ${FILL_FLAGS} \
    -n "${NBINS}" \
    -b "${NSUBINT}" \
    -L "${SUBLENGTH}" \
    --template "${TEMPLATE}" \
    --srcname det_${ARM}_${REP} \
    -f "${FIL}"

rm -f "${RUN_DIR}"/core.*

echo ""
echo "=== task ${IDX} (${ARM} rep ${REP}) complete ==="
echo "End (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"