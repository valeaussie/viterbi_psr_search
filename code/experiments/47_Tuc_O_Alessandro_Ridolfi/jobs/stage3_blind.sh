#!/bin/bash
#SBATCH --job-name=stage3_blind
#SBATCH --output=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/slurm_outputs/slurm-blind-%A_%a.out
#SBATCH --error=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/slurm_outputs/slurm-blind-%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=00:30:00


# ----------------------------------------------------------------------------
# Stage 3 (blind): HMM/Viterbi blind search over a grid of DM trials,
# frequency subbands, and coherent timescales (Nt), lean output.
#
# Grid dimensions (recompute --array if any change):
#   DM:      read from dm_list.txt          (21 trials)
#   Subbands: 100..800 Hz, width 10 Hz, step 9 Hz (1 Hz overlap) (78 subbands)
#   Nt:      read from nt_list.txt          (5 values)
#   Total:   21 * 78 * 5 = 8190 tasks       (array 0-8189)
#
# One Slurm array task = one (DM, subband, Nt) triple.
#
# Index mapping (column-major, DM slowest):
#   sub_idx  = IDX % N_SUB
#   nt_idx   = (IDX // N_SUB) % N_TSFT
#   dm_idx   = IDX // (N_SUB * N_TSFT)
#
# Output layout:
#   stage3_viterbi/blind_v1/DM<XX.XX>/Nt<N>/f0_<F>/
#       blind_Nt<N>_f0_<F>_loglike_curve.dat
#       blind_Nt<N>_f0_<F>_paths.dat
#       blind_Nt<N>_f0_<F>.params
#
# A separate aggregator (aggregate_blind.py) peak-finds across all runs.
# ----------------------------------------------------------------------------

# --- paths (edit if needed) -------------------------------------------------
SETUP="/fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh"
CODE_DIR="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code"
JOBS_DIR="${CODE_DIR}/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs"
VITERBI="${CODE_DIR}/viterbi_pipeline.py"
EXP_DIR="${CODE_DIR}/experiments/47_Tuc_O_Alessandro_Ridolfi/47Tuc_blind_search_v1"

DM_LIST="${JOBS_DIR}/dm_list.txt"
NT_LIST_FILE="${JOBS_DIR}/nt_list.txt"

ROOTNAME="47Tuc"
OUT_BASE="${EXP_DIR}/stage3_viterbi/blind_v1"

# Observation parameters (read from .inf at runtime per DM; kept here as
# fallback reference only -- tsamp and mjd_start are read from the .inf)
TOBS="14395.423083663552"   # seconds; used only to compute Tsft = Tobs/Nt

# Subband grid
BAND_LO=100.0
BAND_HI=800.0
DFB=10.0     # subband width (Hz)
STEP=9.0     # subband step  (Hz; 1 Hz overlap)

# ----------------------------------------------------------------------------
IDX=$(( SLURM_ARRAY_TASK_ID + ${IDX_OFFSET:-0} ))

echo "Start (UTC):  $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Host:         $(hostname)"
echo "Array job:    ${SLURM_ARRAY_JOB_ID:-N/A} task ${IDX}"

# --- read DM and Nt lists ---------------------------------------------------
mapfile -t DM_VALS < <(grep -v '^\s*#' "${DM_LIST}" | grep -v '^\s*$')
mapfile -t NT_VALS < <(grep -v '^\s*#' "${NT_LIST_FILE}" | grep -v '^\s*$')

N_DM=${#DM_VALS[@]}
N_TSFT=${#NT_VALS[@]}

# --- build subband list -----------------------------------------------------
mapfile -t SUBBANDS < <(awk \
    -v lo="${BAND_LO}" -v hi="${BAND_HI}" -v s="${STEP}" \
    'BEGIN{ f=lo; while (f < hi) { printf "%.3f\n", f; f+=s } }')

N_SUB=${#SUBBANDS[@]}
N_TOTAL=$(( N_DM * N_SUB * N_TSFT ))

if [ "${IDX}" -ge "${N_TOTAL}" ]; then
    echo "Array index ${IDX} >= N_TOTAL ${N_TOTAL}; nothing to do."
    exit 0
fi

# --- map array index to (dm_idx, nt_idx, sub_idx) ---------------------------
SUB_IDX=$(( IDX % N_SUB ))
NT_IDX=$(( (IDX / N_SUB) % N_TSFT ))
DM_IDX=$(( IDX / (N_SUB * N_TSFT) ))

DM="${DM_VALS[${DM_IDX}]}"
NT="${NT_VALS[${NT_IDX}]}"
F0_SUB="${SUBBANDS[${SUB_IDX}]}"

DM_TAG=$(printf "DM%05.2f" "${DM}")
TSFT=$(awk "BEGIN { printf \"%.10f\", ${TOBS} / ${NT} }")

echo "  DM           = ${DM}  (${DM_TAG})"
echo "  Nt           = ${NT}  Tsft = ${TSFT} s"
echo "  Subband f0   = ${F0_SUB} Hz"

# --- locate input .dat and .inf for this DM ---------------------------------
BARY_DIR="${EXP_DIR}/stage2_dedisp/${DM_TAG}_bary"
DAT="${BARY_DIR}/${ROOTNAME}_${DM_TAG}_bary.dat"
INF="${BARY_DIR}/${ROOTNAME}_${DM_TAG}_bary.inf"

if [ ! -f "${DAT}" ]; then
    echo "ERROR: .dat not found: ${DAT}" >&2
    exit 1
fi
if [ ! -f "${INF}" ]; then
    echo "ERROR: .inf not found: ${INF}" >&2
    exit 1
fi

# --- read tsamp and mjd_start from .inf -------------------------------------
TSAMP=$(grep -i "width of each time series bin" "${INF}" \
        | awk -F'=' '{gsub(/ /,"",$2); print $2}')
MJD_START=$(grep -i "epoch of observation (mjd)" "${INF}" \
            | awk -F'=' '{gsub(/ /,"",$2); print $2}')

if [ -z "${TSAMP}" ] || [ -z "${MJD_START}" ]; then
    echo "ERROR: could not parse tsamp or mjd_start from ${INF}" >&2
    exit 1
fi
echo "  tsamp        = ${TSAMP} s"
echo "  mjd_start    = ${MJD_START}"

# --- output directory -------------------------------------------------------
F0_FMT=$(printf "%.3f" "${F0_SUB}")
OUT_DIR="${OUT_BASE}/${DM_TAG}/Nt${NT}/f0_${F0_FMT}"
mkdir -p "${OUT_DIR}"
OUT_PREFIX="${OUT_DIR}/blind_Nt${NT}_f0_${F0_FMT}"

# --- environment ------------------------------------------------------------
set +e
source "${SETUP}"
set -e
set -uo pipefail

# --- write params file ------------------------------------------------------
PARAMS="${OUT_PREFIX}.params"
cat > "${PARAMS}" << PARAMEOF
infile     = "${DAT}"
tsamp      = ${TSAMP}
Tsft       = ${TSFT}
Nsft       = -1
f0         = ${F0_SUB}
bw         = ${DFB}
padding_factor = 1
num_harm   = 1
top_paths  = 1
save_delta = False
out_prefix = "${OUT_PREFIX}"
dm         = ${DM}
mjd_start  = ${MJD_START}
PARAMEOF

# --- run Viterbi (lean mode) ------------------------------------------------
echo ""
echo "Running viterbi_pipeline.py --lean ..."
python "${VITERBI}" --params "${PARAMS}" --lean

echo ""
echo "=== task ${IDX} (${DM_TAG} Nt${NT} f0_${F0_FMT}) complete ==="
echo "End (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Output in: ${OUT_DIR}"