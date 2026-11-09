#!/bin/bash
#SBATCH --job-name=viterbi_target
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# ---------------------------------------------------------------------------
# Targeted Viterbi search at a chosen frequency, over a grid of DM trials.
#
# One array task = one DM trial. Everything about the target is passed in
# through the environment, so this one script serves every pulsar and every
# control run without being edited.
#
# Required:
#   TARGET   short name; output goes to redetect_<TARGET>/viterbi_out/
#   F0_HZ    centre frequency of the search band, Hz
#   NT       number of coherent segments
#
# Optional:
#   DMS      space-separated DM trials   (default: 24.30 ... 24.80)
#   BW_HZ    width of the search band    (default: 2.0)
#   LABEL    annotation drawn on the plots (default: TARGET)
#
# Nt is set per submission rather than per array index because memory scales
# as 1/Nt: the segment length is Tobs/Nt, so halving Nt doubles the
# per-segment FFT. One --mem for the whole grid is either wasteful at large
# Nt or fatal at small Nt.
#
# Examples
#   sbatch --array=0-5 --mem=96g --time=00:45:00 \
#     --export=ALL,TARGET=C,F0_HZ=173.708219,NT=2 run_viterbi_target.sh
#
#   sbatch --array=0-5 --mem=12g --time=00:20:00 \
#     --export=ALL,TARGET=E,F0_HZ=282.779107,NT=16,DMS="24.00 24.10 24.20 24.30 24.40 24.50" \
#     run_viterbi_target.sh
#
# Run with --show to print the resolved settings without submitting.
# ---------------------------------------------------------------------------

set -uo pipefail

DMS=${DMS:-"24.30 24.40 24.50 24.60 24.70 24.80"}
BW_HZ=${BW_HZ:-2.0}

SETUP="/fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh"
CODE_DIR="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code"
EXP="${CODE_DIR}/experiments/47_Tuc_O_Alessandro_Ridolfi/47Tuc_blind_search_v1"
SEARCH="${CODE_DIR}/viterbi_search_and_fit.py"
DEDISP="${EXP}/stage2_dedisp"

read -r -a DM_ARR <<< "${DMS}"
N_DM=${#DM_ARR[@]}

if [ "${1:-}" = "--show" ]; then
    echo "TARGET  = ${TARGET:-<unset>}"
    echo "F0_HZ   = ${F0_HZ:-<unset>}"
    echo "NT      = ${NT:-<unset>}"
    echo "BW_HZ   = ${BW_HZ}"
    echo "DMS     = ${DMS}   (${N_DM} tasks)"
    echo "outputs -> ${EXP}/redetect_${TARGET:-<TARGET>}/viterbi_out/"
    echo ""
    echo "submit with --array=0-$(( N_DM - 1 ))"
    exit 0
fi

for v in TARGET F0_HZ NT; do
    if [ -z "${!v:-}" ]; then
        echo "ERROR: ${v} is not set. Pass it with --export=ALL,${v}=..." >&2
        echo "Run '$(basename "$0") --show' to see the resolved settings." >&2
        exit 1
    fi
done

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo "ERROR: no SLURM_ARRAY_TASK_ID. Submit with sbatch --array=..." >&2
    exit 1
fi

LABEL=${LABEL:-${TARGET}}
WORK="${EXP}/redetect_${TARGET}"
OUT_BASE="${WORK}/viterbi_out"

IDX=${SLURM_ARRAY_TASK_ID}
if [ "${IDX}" -ge "${N_DM}" ]; then
    echo "Array index ${IDX} >= ${N_DM}; nothing to do."
    exit 0
fi

DM=${DM_ARR[${IDX}]}
DM_TAG=$(printf "DM%05.2f" "${DM}")

echo "Start (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Host:        $(hostname)"
echo "Memory req:  ${SLURM_MEM_PER_NODE:-unset} MB"
echo "Target:      ${LABEL}  f0 = ${F0_HZ} Hz"
echo "Task ${IDX} of ${N_DM}: ${DM_TAG}  Nt = ${NT}"

DAT="${DEDISP}/${DM_TAG}_bary/47Tuc_${DM_TAG}_bary.dat"
INF="${DEDISP}/${DM_TAG}_bary/47Tuc_${DM_TAG}_bary.inf"

for f in "${DAT}" "${INF}"; do
    if [ ! -f "${f}" ]; then
        echo "ERROR: not found: ${f}" >&2
        ls -d "${DEDISP}"/DM*_bary 2>/dev/null >&2
        exit 1
    fi
done

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
F0_HI=$(awk "BEGIN { printf \"%.6f\", ${F0_HZ} + ${BW_HZ}/2.0 }")
BINW=$(awk "BEGIN { printf \"%.3e\", 1.0/${TSFT} }")

echo "  tsamp       = ${TSAMP} s"
echo "  mjd_start   = ${MJD_START}   (barycentric; the filterbank header differs)"
echo "  nsamp       = ${NSAMP}   (prepsubband pads; trimmed at read time)"
echo "  Tobs        = ${TOBS} s"
echo "  Tcoh        = ${TSFT} s   (Nt = ${NT})"
echo "  bin width   = ${BINW} Hz"
echo "  samples/seg = ${NPSEG}   <- this drives memory"
echo "  band        = ${F0_LO} to ${F0_HI} Hz"

set +e
source "${SETUP}"
set -e
set -uo pipefail

OUT_DIR="${OUT_BASE}/${DM_TAG}_Nt${NT}"
mkdir -p "${OUT_DIR}"
PREFIX="${OUT_DIR}/redetect_${TARGET}_${DM_TAG}_Nt${NT}"

cd "${OUT_DIR}"
python3 "${SEARCH}" \
    --infile "${DAT}" \
    --tsamp "${TSAMP}" \
    --Tsft "${TSFT}" --Nsft "${NT}" \
    --f0 "${F0_LO}" --bw "${BW_HZ}" \
    --poly-order 5 --plot-path \
    --dm "${DM}" --mjd-start "${MJD_START}" \
    --known-label "${LABEL}" \
    --out_prefix "${PREFIX}" \
    2>&1 | tee "${PREFIX}.log"

RC=${PIPESTATUS[0]}
echo "=== ${TARGET} task ${IDX} (${DM_TAG} Nt${NT}) exit ${RC} ==="
echo "End (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"
ls -l "${OUT_DIR}"
exit "${RC}"