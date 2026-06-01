#!/bin/bash
#SBATCH --job-name=stage3b_fit
#SBATCH --output=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/slurm_outputs/slurm-fit-%A_%a.out
#SBATCH --error=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/slurm_outputs/slurm-fit-%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g
#SBATCH --time=00:30:00
#SBATCH --array=0-376   # 377 candidates (0-indexed); adjust if dedup count changes

# ----------------------------------------------------------------------------
# Stage 3b (fit): run viterbi_search_and_fit.py on each deduplicated candidate
# from candidates_dedup.csv.
#
# Each task:
#   1. Calls parse_dedup_row.py to safely parse the candidate's fields from
#      candidates_dedup.csv (handles the nt_values multi-column quirk).
#   2. Locates the existing .params file from the blind run for that
#      (Nt_best, subband_f0) pair.
#   3. Writes a new .params file with a 10 Hz window centred on peak_freq_hz,
#      Tsft from Nt_best, and output prefix pointing to the candidate directory.
#   4. Runs viterbi_search_and_fit.py --params --plot-path.
#
# Output layout
# -------------
#   stage3_viterbi/candidates/
#     cand_000/
#       cand_000.params
#       cand_000_psrfold.candfile
#       cand_000_psrfold.candfile.info
#       cand_000_fit_summary.dat
#       cand_000_spectrogram.png
#       cand_000_track.dat
#       cand_000_loglikes.png
#     cand_001/
#       ...
#
# Prerequisites
# -------------
#   - Stage 3 blind search complete (stage3_blind.sh).
#   - aggregate_blind.py run, producing candidates_dedup.csv.
#   - parse_dedup_row.py present in CODE_DIR (alongside this script).
# ----------------------------------------------------------------------------

# --- paths (edit if needed) -------------------------------------------------
SETUP="/fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh"
CODE_DIR="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code"
FIT_SCRIPT="${CODE_DIR}/viterbi_search_and_fit.py"
PARSE_SCRIPT="${CODE_DIR}/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/parse_dedup_row.py"
EXP_DIR="${CODE_DIR}/experiments/47_Tuc_O_Alessandro_Ridolfi/47Tuc_blind_search_v1"

DEDUP_CSV="${EXP_DIR}/stage3_viterbi/blind_v1/candidates_dedup.csv"
BLIND_BASE="${EXP_DIR}/stage3_viterbi/blind_v1"
CAND_DIR="${EXP_DIR}/stage3_viterbi/candidates"

# Observation total length in seconds (must match stage3_blind.sh)
TOBS="14395.423083663552"

# Candidate window width in Hz (centred on peak_freq_hz)
BW_CAND="10.0"

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

# --- parse candidate fields safely via Python helper ------------------------
PARSED=$(python "${PARSE_SCRIPT}" --csv "${DEDUP_CSV}" --idx "${IDX}")
if [ $? -ne 0 ]; then
    echo "ERROR: parse_dedup_row.py failed for index ${IDX}"
    exit 1
fi

# Evaluate into shell variables: PEAK_FREQ, PEAK_LL, MULT, NT_VALS,
#                                NT_BEST, SUB_F0, KNOWN
eval "${PARSED}"

echo "  peak_freq    = ${PEAK_FREQ} Hz"
echo "  peak_loglike = ${PEAK_LL}"
echo "  Nt_best      = ${NT_BEST}"
echo "  subband_f0   = ${SUB_F0} Hz"
echo "  multiplicity = ${MULT}  nt_values = ${NT_VALS}"
echo "  known_match  = ${KNOWN:-none}"

# --- locate the blind-run .params file for (Nt_best, subband_f0) ------------
SUB_F0_FMT=$(printf "%.3f" "${SUB_F0}")
BLIND_RUN_DIR="${BLIND_BASE}/Nt${NT_BEST}/f0_${SUB_F0_FMT}"
BLIND_PARAMS="${BLIND_RUN_DIR}/blind_Nt${NT_BEST}_f0_${SUB_F0_FMT}.params"

if [ ! -f "${BLIND_PARAMS}" ]; then
    echo "ERROR: blind params file not found: ${BLIND_PARAMS}"
    exit 1
fi

# --- compute Tsft = Tobs / Nt_best ------------------------------------------
TSFT=$(awk "BEGIN { printf \"%.10f\", ${TOBS} / ${NT_BEST} }")
echo "  Tsft         = ${TSFT} s  (Tobs / Nt = ${TOBS} / ${NT_BEST})"

# --- candidate window: centre BW_CAND on peak_freq --------------------------
F0_CAND=$(awk "BEGIN { printf \"%.10f\", ${PEAK_FREQ} - ${BW_CAND}/2.0 }")
F1_CAND=$(awk "BEGIN { printf \"%.6f\",  ${PEAK_FREQ} + ${BW_CAND}/2.0 }")
echo "  window       = [${F0_CAND}, ${F1_CAND}] Hz"

# --- extract global obs params from the blind params file -------------------
INFILE=$(   grep -m1 '^infile'     "${BLIND_PARAMS}" | awk -F'=' '{gsub(/[ "]/,"",$2); print $2}')
TSAMP=$(    grep -m1 '^tsamp'      "${BLIND_PARAMS}" | awk -F'=' '{gsub(/ /,"",$2);   print $2}')
DM=$(       grep -m1 '^dm'         "${BLIND_PARAMS}" | awk -F'=' '{gsub(/ /,"",$2);   print $2}')
MJD_START=$(grep -m1 '^mjd_start'  "${BLIND_PARAMS}" | awk -F'=' '{gsub(/ /,"",$2);   print $2}')

echo "  infile       = ${INFILE}"
echo "  tsamp        = ${TSAMP} s"
echo "  dm           = ${DM}"
echo "  mjd_start    = ${MJD_START}"

# --- create output directory ------------------------------------------------
OUT_DIR="${CAND_DIR}/${CAND_ID}"
mkdir -p "${OUT_DIR}"
OUT_PREFIX="${OUT_DIR}/${CAND_ID}"

# --- write per-candidate params file ----------------------------------------
PARAMS="${OUT_DIR}/${CAND_ID}.params"

cat > "${PARAMS}" << PARAMEOF
infile = "${INFILE}"
tsamp = ${TSAMP}
Tsft = ${TSFT}
Nsft = -1
f0 = ${F0_CAND}
bw = ${BW_CAND}
padding_factor = 1
out_prefix = "${OUT_PREFIX}"
plot_path = True
top_paths = 1
save_delta = False
num_harm = 1
dm = ${DM}
mjd_start = ${MJD_START}
PARAMEOF

echo ""
echo "Params written: ${PARAMS}"
echo "Running viterbi_search_and_fit.py ..."
echo ""

# --- run the fit ------------------------------------------------------------
# --spec-flo/fhi zooms the spectrogram plot tightly around the candidate freq
SPEC_LO=$(awk "BEGIN { printf \"%.6f\", ${PEAK_FREQ} - 0.05 }")
SPEC_HI=$(awk "BEGIN { printf \"%.6f\", ${PEAK_FREQ} + 0.05 }")

python "${FIT_SCRIPT}" \
    --params  "${PARAMS}" \
    --plot-path \
    --spec-flo "${SPEC_LO}" \
    --spec-fhi "${SPEC_HI}"

echo ""
echo "=== task ${IDX} (${CAND_ID}) complete ==="
echo "End (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Output in: ${OUT_DIR}"