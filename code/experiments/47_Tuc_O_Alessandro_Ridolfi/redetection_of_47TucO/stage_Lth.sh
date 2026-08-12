#!/bin/bash
#SBATCH --job-name=Lth_grid
#SBATCH --output=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/redetection_of_47TucO/slurm_outputs/slurm-Lth-%A_%a.out
#SBATCH --error=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/redetection_of_47TucO/slurm_outputs/slurm-Lth-%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=01:00:00
#SBATCH --array=0-14

# ----------------------------------------------------------------------------
# Detection threshold L_th(N_T) for the coherent-timescale grid of 47 Tuc O.
# One array task per N_T. Each task:
#   1. Reads N_Q (= data rows of that run's *_loglikes.dat).
#   2. Runs compute_Lth.py (exponential-tail method, O'Leary et al. 2026),
#      writing one CSV row to Lth_grid_parts/Lth_Nt<NT>.csv.
# After the array finishes, run collate_Lth.py to build Lth_grid.csv.
#
# Before submitting:  mkdir -p slurm_outputs Lth_grid_parts
# ----------------------------------------------------------------------------

# --- paths (edit if needed) -------------------------------------------------
SETUP="/fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh"
REDET_DIR="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/redetection_of_47TucO"
COMPUTE="${REDET_DIR}/compute_Lth.py"
PARTS_DIR="${REDET_DIR}/Lth_grid_parts"
ROOT="47TucO_redetect"

# --- fixed analysis choices -------------------------------------------------
NREAL=10000
ALPHA=0.1
NBOOT=200

# --- grid (must match the search grid) --------------------------------------
NTS=(16 20 24 28 30 31 32 40 48 56 64 128 256 512 1024)
TSFTS=(899.713885 719.771093 599.809282 514.122187 479.847395 464.368420 \
       449.856943 359.885508 299.904603 257.061055 224.928433 112.464217 \
       56.232108 28.116016 14.058008)

# ----------------------------------------------------------------------------
IDX=${SLURM_ARRAY_TASK_ID}
NT=${NTS[$IDX]}
TC=${TSFTS[$IDX]}
LLFILE="${REDET_DIR}/${ROOT}_Nt${NT}_loglikes.dat"
OUT="${PARTS_DIR}/Lth_Nt${NT}.csv"

echo "Start (UTC):  $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Host:         $(hostname)"
echo "Array job:    ${SLURM_ARRAY_JOB_ID:-N/A} task ${IDX}"
echo "N_T:          ${NT}   T_coh: ${TC} s"

# --- environment ------------------------------------------------------------
set +e
source "${SETUP}"
set -e
set -uo pipefail

mkdir -p "${PARTS_DIR}"

if [ ! -f "${LLFILE}" ]; then
    echo "ERROR: loglikes file not found: ${LLFILE}" >&2
    exit 1
fi

NBIN=$(( $(wc -l < "${LLFILE}") - 1 ))
echo "N_Q (Nbin):   ${NBIN}"

rm -f "${OUT}"   # start clean so a re-run does not append a duplicate row

python "${COMPUTE}" \
    --NT "${NT}" \
    --Nbin "${NBIN}" \
    --Nreal "${NREAL}" \
    --alpha-prime "${ALPHA}" \
    --nboot "${NBOOT}" \
    --Tcoh "${TC}" \
    --seed "${NT}" \
    --out "${OUT}"

echo "Wrote ${OUT}"
echo "End (UTC):    $(date -u +%Y-%m-%dT%H:%M:%SZ)"