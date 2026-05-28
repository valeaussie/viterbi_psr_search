#!/bin/bash
#SBATCH --job-name=stage3_blind
#SBATCH --output=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/slurm_outputs/slurm-blind-%A_%a.out
#SBATCH --error=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/slurm_outputs/slurm-blind-%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=00:30:00
#SBATCH --array=0-389

# ----------------------------------------------------------------------------
# Stage 3 (blind): HMM/Viterbi blind search over a grid of frequency subbands
# and coherent timescales (Tsft), single DM, lean output.
#
# Grid:
#   - Subbands: 100..800 Hz, width 10 Hz, step 9 Hz (1 Hz overlap).
#   - Tsft (Nt) grid: Nt in {16,32,64,128,256}, Tsft = Tobs/Nt.
#
# One Slurm array task = one (subband, Tsft) pair. Each task runs
# viterbi_pipeline.py --lean, which saves only the loglike-vs-frequency curve
# and the top paths (no fits, no candfile, no plots). A separate aggregator
# (aggregate_blind.py) peak-finds across all runs afterwards.
#
# Array size must equal n_subbands * n_Tsft - 1. With 78 subbands x 5 Tsft
# = 390 runs, use --array=0-389. (Recompute if you change the grid.)
# ----------------------------------------------------------------------------

# --- paths (edit if needed) -------------------------------------------------
SETUP="/fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh"
CODE_DIR="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code"
VITERBI="${CODE_DIR}/viterbi_pipeline.py"
EXP_DIR="${CODE_DIR}/experiments/47_Tuc_O_Alessandro_Ridolfi/47Tuc_blind_search_v1"

INFILE="${CODE_DIR}/experiments/47_Tuc_O_Alessandro_Ridolfi/data/J0000-00_cfbf00000_Plan1_1_DM24.36.dat"
TSAMP="7.65607476635514e-05"
DM="24.356441"
MJD_START="59391.125996527684038"

OUT_BASE="${EXP_DIR}/stage3_viterbi/blind_v1"

# --- grid definition (must match the array size above) ----------------------
# Subbands
BAND_LO=100.0
BAND_HI=800.0
DFB=10.0
STEP=9.0

# Tsft grid (exact: Tobs/Nt with Tobs = 14395.423083663552 s)
# Nt:    16        32        64        128       256
TSFT_LIST=(899.7139 449.8570 224.9285 112.4642 56.2321)
NT_LIST=(16 32 64 128 256)

# --- build the full subband edge list (bash float loop via awk) -------------
mapfile -t SUBBANDS < <(awk -v lo="${BAND_LO}" -v hi="${BAND_HI}" -v w="${DFB}" -v s="${STEP}" \
    'BEGIN{ f=lo; while (f < hi) { printf "%.3f\n", f; f+=s } }')

N_SUB=${#SUBBANDS[@]}
N_TSFT=${#TSFT_LIST[@]}
N_TOTAL=$(( N_SUB * N_TSFT ))

# --- map this array index to (subband, Tsft) -------------------------------
IDX=${SLURM_ARRAY_TASK_ID}
if [ "${IDX}" -ge "${N_TOTAL}" ]; then
    echo "Array index ${IDX} >= N_TOTAL ${N_TOTAL}; nothing to do."
    exit 0
fi

SUB_I=$(( IDX % N_SUB ))
TSFT_I=$(( IDX / N_SUB ))

F0=${SUBBANDS[${SUB_I}]}
TSFT=${TSFT_LIST[${TSFT_I}]}
NT=${NT_LIST[${TSFT_I}]}

# --- header -----------------------------------------------------------------
echo "Start (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Host:        $(hostname)"
echo "Array job:   ${SLURM_ARRAY_JOB_ID:-N/A} task ${IDX}"
echo "Grid:        ${N_SUB} subbands x ${N_TSFT} Tsft = ${N_TOTAL} runs"
echo "This task:   f0=${F0} Hz, bw=${DFB} Hz, Tsft=${TSFT} s (Nt=${NT})"

# --- environment ------------------------------------------------------------
set +e
source "${SETUP}"
set -e
set -uo pipefail

# --- per-run output dir and params file -------------------------------------
RUN_DIR="${OUT_BASE}/Nt${NT}/f0_${F0}"
mkdir -p "${RUN_DIR}"
cd "${RUN_DIR}"

OUT_PREFIX="blind_Nt${NT}_f0_${F0}"
PARAMS="${RUN_DIR}/${OUT_PREFIX}.params"

cat > "${PARAMS}" << PARAMEOF
infile = "${INFILE}"
tsamp = ${TSAMP}
Tsft = ${TSFT}
Nsft = -1
f0 = ${F0}
bw = ${DFB}
padding_factor = 1
out_prefix = "${OUT_PREFIX}"
plot_path = False
top_paths = 1
save_delta = False
num_harm = 1
dm = ${DM}
mjd_start = ${MJD_START}
PARAMEOF

# --- run the lean Viterbi search --------------------------------------------
echo ""
echo "Running viterbi_pipeline.py --lean ..."
python "${VITERBI}" --params "${PARAMS}" --lean

echo ""
echo "=== task ${IDX} complete (f0=${F0}, Nt=${NT}) ==="
echo "End (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Loglike curve: ${RUN_DIR}/${OUT_PREFIX}_loglike_curve.dat"