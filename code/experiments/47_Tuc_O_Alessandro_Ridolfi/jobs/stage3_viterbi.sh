#!/bin/bash
#SBATCH --job-name=stage3_viterbi
#SBATCH --output=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/slurm_outputs/slurm-viterbi-%j.out
#SBATCH --error=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/slurm_outputs/slurm-viterbi-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8g
#SBATCH --time=00:30:00

MODE="${1:-}"
if [ "${MODE}" != "bary" ] && [ "${MODE}" != "topo" ]; then
    echo "ERROR: first argument must be 'bary' or 'topo'." >&2
    echo "Usage: sbatch stage3_viterbi.sh [bary|topo]" >&2
    exit 1
fi

SETUP="/fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh"
CODE_DIR="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code"
VITERBI="${CODE_DIR}/viterbi_pipeline.py"
GENPARAMS="${CODE_DIR}/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/gen_search_params.py"

EXP_DIR="${CODE_DIR}/experiments/47_Tuc_O_Alessandro_Ridolfi/47Tuc_blind_search_v1"
DATA_DIR="${CODE_DIR}/experiments/47_Tuc_O_Alessandro_Ridolfi/data"

DM="24.356441"
DM_TAG="DM24.356441"

STAGE2_DIR="${EXP_DIR}/stage2_dedisp/${DM_TAG}_${MODE}"
DAT="${STAGE2_DIR}/47Tuc_${DM_TAG}_${MODE}.dat"
INF="${STAGE2_DIR}/47Tuc_${DM_TAG}_${MODE}.inf"

PARAMS="${DATA_DIR}/search_${MODE}.params"
OUT_PREFIX="exp_01_${MODE}"
OUT_DIR="${EXP_DIR}/stage3_viterbi/${DM_TAG}_${MODE}"

TSFT=450

echo "Start (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Host:        $(hostname)"
echo "Job ID:      ${SLURM_JOB_ID:-N/A}"
echo "Mode:        ${MODE}"
echo "Input .dat:  ${DAT}"
echo "Tsft:        ${TSFT}"

set +e
source "${SETUP}"
set -e
set -uo pipefail

[ -f "${DAT}" ] || { echo "ERROR: .dat not found: ${DAT}" >&2; exit 1; }
[ -f "${INF}" ] || { echo "ERROR: .inf not found: ${INF}" >&2; exit 1; }
[ -f "${VITERBI}" ] || { echo "ERROR: viterbi code not found: ${VITERBI}" >&2; exit 1; }
[ -f "${GENPARAMS}" ] || { echo "ERROR: gen_search_params.py not found: ${GENPARAMS}" >&2; exit 1; }

mkdir -p "${OUT_DIR}"

echo ""
echo "Generating params file: ${PARAMS}"
python "${GENPARAMS}" \
    --inf "${INF}" \
    --dat "${DAT}" \
    --out "${PARAMS}" \
    --tsft "${TSFT}" \
    --out-prefix "${OUT_PREFIX}" \
    --dm "${DM}"

cd "${OUT_DIR}"

echo ""
echo "Running Viterbi search (mode=${MODE}, Tsft=${TSFT})..."
echo ""
python "${VITERBI}" --params "${PARAMS}"

echo ""
echo "=== Stage 3 (${MODE}, Tsft=${TSFT}) complete ==="
echo "End (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Outputs in: ${OUT_DIR}"
ls -lh "${OUT_DIR}"