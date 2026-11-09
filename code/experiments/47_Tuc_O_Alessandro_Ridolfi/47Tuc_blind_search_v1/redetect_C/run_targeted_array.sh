#!/bin/bash
#SBATCH --job-name=redetect
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00

# ---------------------------------------------------------------------------
# Targeted Viterbi search for one known pulsar.
#
# One array task = one DM trial. Nt is fixed per submission and passed in
# through the environment, because memory scales as 1/Nt: the segment length
# is Tobs/Nt, so halving Nt doubles the per-segment FFT. A single --mem for
# the whole grid is either wasteful at large Nt or fatal at small Nt.
#
# Submit one array per Nt, with memory to match:
#
#   sbatch --array=0-5 --mem=16g  --export=ALL,NT=32 run_targeted_array.sh
#   sbatch --array=0-5 --mem=32g  --export=ALL,NT=16 run_targeted_array.sh
#   sbatch --array=0-5 --mem=64g  --export=ALL,NT=8  run_targeted_array.sh
#   sbatch --array=0-5 --mem=128g --export=ALL,NT=4  run_targeted_array.sh
#
# Those figures are starting guesses, not measurements. Once one task has
# run, read the real requirement and trim them, since an oversized request
# delays scheduling:
#   sacct -j <jobid> --format=JobID,Elapsed,MaxRSS,State
#
# Run with --show-grid to see the submission lines without submitting.
# ---------------------------------------------------------------------------

set -uo pipefail

# --- target -----------------------------------------------------------------
PULSAR="C"
KNOWN_LABEL="J0024-7204C"
F0_HZ=173.708219          # catalogue spin frequency, Hz
BW_HZ=2.0                 # width of the search window, centred on F0_HZ

# --- DM trials (one array task each) ----------------------------------------
DMS=(24.30 24.40 24.50 24.60 24.70 24.80)

# --- paths ------------------------------------------------------------------
SETUP="/fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh"
CODE_DIR="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code"
EXP="${CODE_DIR}/experiments/47_Tuc_O_Alessandro_Ridolfi/47Tuc_blind_search_v1"
SEARCH="${CODE_DIR}/viterbi_search_and_fit.py"
DEDISP="${EXP}/stage2_dedisp"
WORK="${EXP}/redetect_${PULSAR}"
OUT_BASE="${WORK}/viterbi_out"

N_DM=${#DMS[@]}

if [ "${1:-}" = "--show-grid" ]; then
    echo "Target : ${KNOWN_LABEL}  f0 = ${F0_HZ} Hz"
    echo "DM     : ${DMS[*]}   (${N_DM} tasks)"
    echo ""
    echo "Submit one array per Nt, memory scaled to the segment length:"
    echo "  sbatch --array=0-$(( N_DM - 1 )) --mem=16g  --export=ALL,NT=32 $(basename "$0")"
    echo "  sbatch --array=0-$(( N_DM - 1 )) --mem=32g  --export=ALL,NT=16 $(basename "$0")"
    echo "  sbatch --array=0-$(( N_DM - 1 )) --mem=64g  --export=ALL,NT=8  $(basename "$0")"
    echo "  sbatch --array=0-$(( N_DM - 1 )) --mem=128g --export=ALL,NT=4  $(basename "$0")"
    exit 0
fi

if [ -z "${NT:-}" ]; then
    echo "ERROR: NT is not set. Pass it in with --export=ALL,NT=<value>." >&2
    echo "Run '$(basename "$0") --show-grid' for the submission lines." >&2
    exit 1
fi

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo "ERROR: no SLURM_ARRAY_TASK_ID. Submit with sbatch --array=..." >&2
    exit 1
fi

IDX=${SLURM_ARRAY_TASK_ID}
if [ "${IDX}" -ge "${N_DM}" ]; then
    echo "Array index ${IDX} >= ${N_DM}; nothing to do."
    exit 0
fi

DM=${DMS[${IDX}]}
DM_TAG=$(printf "DM%05.2f" "${DM}")

echo "Start (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Host:        $(hostname)"
echo "Memory req:  ${SLURM_MEM_PER_NODE:-unset} MB"
echo "Target:      ${KNOWN_LABEL}  f0 = ${F0_HZ} Hz"
echo "Task ${IDX}: ${DM_TAG}  Nt = ${NT}"

DAT="${DEDISP}/${DM_TAG}_bary/47Tuc_${DM_TAG}_bary.dat"
INF="${DEDISP}/${DM_TAG}_bary/47Tuc_${DM_TAG}_bary.inf"

for f in "${DAT}" "${INF}"; do
    if [ ! -f "${f}" ]; then
        echo "ERROR: not found: ${f}" >&2
        ls -d "${DEDISP}"/DM*_bary 2>/dev/null >&2
        exit 1
    fi
done

# --- observation parameters, read from the .inf -----------------------------
TSAMP=$(grep -i "width of each time series bin" "${INF}" \
        | awk -F'=' '{gsub(/ /,"",$2); print $2}')
MJD_START=$(grep -i "epoch of observation (mjd)" "${INF}" \
            | awk -F'=' '{gsub(/ /,"",$2); print $2}')
NSAMP=$(grep -i "number of bins in the time series" "${INF}" \
        | awk -F'=' '{gsub(/ /,"",$2); print $2}')

if [ -z "${TSAMP}" ] || [ -z "${MJD_START}" ] || [ -z "${NSAMP}" ]; then
    echo "ERROR: could not parse tsamp, mjd_start or nsamp from ${INF}" >&2
    exit 1
fi

TOBS=$(awk "BEGIN { printf \"%.10f\", ${NSAMP} * ${TSAMP} }")
TSFT=$(awk "BEGIN { printf \"%.10f\", ${TOBS} / ${NT} }")
NPSEG=$(awk "BEGIN { printf \"%.0f\", ${NSAMP} / ${NT} }")
F0_LO=$(awk "BEGIN { printf \"%.6f\", ${F0_HZ} - ${BW_HZ}/2.0 }")

echo "  tsamp       = ${TSAMP} s"
echo "  mjd_start   = ${MJD_START}"
echo "  nsamp       = ${NSAMP}   (prepsubband pads; the filterbank is shorter)"
echo "  Tobs        = ${TOBS} s"
echo "  Tcoh        = ${TSFT} s   (Nt = ${NT})"
echo "  samples/seg = ${NPSEG}   <- this is what drives memory"
echo "  band        = ${F0_LO} to $(awk "BEGIN{printf \"%.6f\", ${F0_LO}+${BW_HZ}}") Hz"

# --- environment ------------------------------------------------------------
set +e
source "${SETUP}"
set -e
set -uo pipefail

OUT_DIR="${OUT_BASE}/${DM_TAG}_Nt${NT}"
mkdir -p "${OUT_DIR}"
PREFIX="${OUT_DIR}/redetect_${PULSAR}_${DM_TAG}_Nt${NT}"

cd "${OUT_DIR}"
python3 "${SEARCH}" \
    --infile "${DAT}" \
    --tsamp "${TSAMP}" \
    --Tsft "${TSFT}" --Nsft "${NT}" \
    --f0 "${F0_LO}" --bw "${BW_HZ}" \
    --poly-order 5 --plot-path \
    --dm "${DM}" --mjd-start "${MJD_START}" \
    --known-label "${KNOWN_LABEL}" \
    --out_prefix "${PREFIX}" \
    2>&1 | tee "${PREFIX}.log"

RC=${PIPESTATUS[0]}
echo "=== task ${IDX} (${DM_TAG} Nt${NT}) exit ${RC} ==="
echo "End (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"
ls -l "${OUT_DIR}"
exit "${RC}"