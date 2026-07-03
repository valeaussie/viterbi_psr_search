#!/bin/bash
#SBATCH --job-name=stage3b_fit
#SBATCH --output=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/slurm_outputs/slurm-fit-%A_%a.out
#SBATCH --error=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/slurm_outputs/slurm-fit-%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=00:10:00
#SBATCH --array=0-199   # adjust upper bound to len(candidates_dedup.csv) - 2

# ----------------------------------------------------------------------------
# Stage 3b (fit): run viterbi_search_and_fit.py on each deduplicated candidate
# from candidates_dedup.csv.
#
# Each task:
#   1. Calls parse_dedup_row.py to read all candidate fields, including
#      DM_BEST (the DM at which the peak loglike was highest).
#   2. Locates the dedispersed .dat and .inf for DM_BEST from stage2_dedisp.
#   3. Reads tsamp and mjd_start from the .inf (never hardcoded).
#   4. Locates the blind-run .params file for (DM_BEST, Nt_best, subband_f0).
#   5. Writes a per-candidate params file with a 10 Hz window centred on
#      peak_freq_hz, Tsft = Tobs/Nt_best, and dm = DM_BEST.
#   6. Runs viterbi_search_and_fit.py.
#
# Output layout:
#   stage3_viterbi/candidates/
#     cand_000/
#       cand_000.params
#       cand_000_psrfold.candfile
#       cand_000_psrfold.candfile.info
#       cand_000_fit_summary.dat
#       cand_000_spectrogram.png
#       cand_000_track.dat
#       cand_000_loglikes.png
#
# Prerequisites:
#   Stage 3 blind complete, aggregate_blind.py run.
# ----------------------------------------------------------------------------

# --- paths (edit if needed) -------------------------------------------------
SETUP="/fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh"
CODE_DIR="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code"
JOBS_DIR="${CODE_DIR}/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs"
FIT_SCRIPT="${CODE_DIR}/viterbi_search_and_fit.py"
PARSE_SCRIPT="${JOBS_DIR}/parse_dedup_row.py"
EXP_DIR="${CODE_DIR}/experiments/47_Tuc_O_Alessandro_Ridolfi/47Tuc_blind_search_v1"

DEDUP_CSV="${EXP_DIR}/stage3_viterbi/blind_v1/candidates_dedup.csv"
BLIND_BASE="${EXP_DIR}/stage3_viterbi/blind_v1"
CAND_DIR="${EXP_DIR}/stage3_viterbi/candidates"

ROOTNAME="47Tuc"
TOBS="14395.423083663552"   # seconds; used to compute Tsft = Tobs / Nt_best
BW_CAND="2.0"              # candidate window width in Hz

# ----------------------------------------------------------------------------
IDX=${SLURM_ARRAY_TASK_ID}
CAND_ID=$(printf "cand_%03d" "${IDX}")

echo "Start (UTC):  $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Host:         $(hostname)"
echo "Array job:    ${SLURM_ARRAY_JOB_ID:-N/A} task ${IDX}"
echo "Candidate:    ${CAND_ID}"

# --- environment ------------------------------------------------------------
set +e
source "${SETUP}"
set -e
set -uo pipefail

# --- parse candidate fields -------------------------------------------------
PARSED=$(python "${PARSE_SCRIPT}" --csv "${DEDUP_CSV}" --idx "${IDX}")
if [ $? -ne 0 ]; then
    echo "ERROR: parse_dedup_row.py failed for index ${IDX}"
    exit 1
fi

PEAK_FREQ=$(echo "${PARSED}" | grep '^PEAK_FREQ=' | cut -d= -f2)
PEAK_LL=$(echo "${PARSED}"   | grep '^PEAK_LL='   | cut -d= -f2)
DM_BEST=$(echo "${PARSED}"   | grep '^DM_BEST='   | cut -d= -f2)
DM_COUNT=$(echo "${PARSED}"  | grep '^DM_COUNT='  | cut -d= -f2)
DM_VALUES=$(echo "${PARSED}" | grep '^DM_VALUES=' | cut -d= -f2)
MULT=$(echo "${PARSED}"      | grep '^MULT='      | cut -d= -f2)
NT_VALS=$(echo "${PARSED}"   | grep '^NT_VALS='   | cut -d= -f2)
NT_BEST=$(echo "${PARSED}"   | grep '^NT_BEST='   | cut -d= -f2)
SUB_F0=$(echo "${PARSED}"    | grep '^SUB_F0='    | cut -d= -f2)
KNOWN=$(echo "${PARSED}"     | grep '^KNOWN='     | cut -d= -f2)
# Now available: PEAK_FREQ, PEAK_LL, DM_BEST, DM_COUNT, DM_VALUES,
#                MULT, NT_VALS, NT_BEST, SUB_F0, KNOWN

echo "  peak_freq    = ${PEAK_FREQ} Hz"
echo "  peak_loglike = ${PEAK_LL}"
echo "  dm_best      = ${DM_BEST}  (recovered at ${DM_COUNT} DM trials: ${DM_VALUES})"
echo "  Nt_best      = ${NT_BEST}"
echo "  subband_f0   = ${SUB_F0} Hz"
echo "  multiplicity = ${MULT}  nt_values = ${NT_VALS}"
echo "  known_match  = ${KNOWN:-none}"

# --- locate .dat and .inf for DM_BEST ---------------------------------------
DM_TAG=$(printf "DM%05.2f" "${DM_BEST}")
BARY_DIR="${EXP_DIR}/stage2_dedisp/${DM_TAG}_bary"
DAT="${BARY_DIR}/${ROOTNAME}_${DM_TAG}_bary.dat"
INF="${BARY_DIR}/${ROOTNAME}_${DM_TAG}_bary.inf"

if [ ! -f "${DAT}" ]; then
    echo "ERROR: .dat not found for DM_BEST=${DM_BEST}: ${DAT}" >&2
    exit 1
fi
if [ ! -f "${INF}" ]; then
    echo "ERROR: .inf not found for DM_BEST=${DM_BEST}: ${INF}" >&2
    exit 1
fi

# --- read tsamp and mjd_start from .inf (never hardcoded) -------------------
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

# --- locate blind params file for (DM_BEST, Nt_best, subband_f0) ------------
SUB_F0_FMT=$(printf "%.3f" "${SUB_F0}")
BLIND_RUN_DIR="${BLIND_BASE}/${DM_TAG}/Nt${NT_BEST}/f0_${SUB_F0_FMT}"
BLIND_PARAMS="${BLIND_RUN_DIR}/blind_Nt${NT_BEST}_f0_${SUB_F0_FMT}.params"

if [ ! -f "${BLIND_PARAMS}" ]; then
    echo "ERROR: blind params not found: ${BLIND_PARAMS}" >&2
    exit 1
fi

# --- compute Tsft -----------------------------------------------------------
TSFT=$(awk "BEGIN { printf \"%.10f\", ${TOBS} / ${NT_BEST} }")
echo "  Tsft         = ${TSFT} s  (Tobs / Nt = ${TOBS} / ${NT_BEST})"

# --- candidate frequency window ---------------------------------------------
F0_CAND=$(awk "BEGIN { printf \"%.10f\", ${PEAK_FREQ} - ${BW_CAND}/2.0 }")
echo "  window       = [${F0_CAND}, +${BW_CAND}] Hz"

# --- output directory -------------------------------------------------------
OUT_DIR="${CAND_DIR}/${CAND_ID}"
mkdir -p "${OUT_DIR}"
OUT_PREFIX="${OUT_DIR}/${CAND_ID}"

# --- write per-candidate params file ----------------------------------------
PARAMS="${OUT_DIR}/${CAND_ID}.params"
cat > "${PARAMS}" << PARAMEOF
infile     = "${DAT}"
tsamp      = ${TSAMP}
Tsft       = ${TSFT}
Nsft       = -1
f0         = ${F0_CAND}
bw         = ${BW_CAND}
padding_factor = 1
num_harm   = 1
top_paths  = 1
save_delta = False
out_prefix = "${OUT_PREFIX}"
plot_path  = True
dm         = ${DM_BEST}
mjd_start  = ${MJD_START}
PARAMEOF

echo ""
echo "Params written: ${PARAMS}"
echo "Running viterbi_search_and_fit.py ..."
echo ""

# --- run the fit ------------------------------------------------------------
SPEC_LO=$(awk "BEGIN { printf \"%.6f\", ${PEAK_FREQ} - 0.05 }")
SPEC_HI=$(awk "BEGIN { printf \"%.6f\", ${PEAK_FREQ} + 0.05 }")

python "${FIT_SCRIPT}" \
    --params   "${PARAMS}" \
    --plot-path \
    --spec-flo "${SPEC_LO}" \
    --spec-fhi "${SPEC_HI}"

echo ""
echo "=== task ${IDX} (${CAND_ID}) complete ==="
echo "End (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Output in: ${OUT_DIR}"