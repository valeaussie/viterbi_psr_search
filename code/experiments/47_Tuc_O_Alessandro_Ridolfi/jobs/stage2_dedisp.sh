#!/bin/bash
#SBATCH --job-name=stage2_dedisp
#SBATCH --output=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/slurm_outputs/slurm-dedisp-%j.out
#SBATCH --error=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/slurm_outputs/slurm-dedisp-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16g
#SBATCH --tmp=40g
#SBATCH --time=02:00:00

# ----------------------------------------------------------------------------
# Stage 2: incoherent dedispersion + barycentring with PRESTO prepdata.
#
# Produces TWO time series from the cleaned filterbank, at a single DM:
#   - barycentred  (prepdata default)
#   - topocentric  (prepdata -nobary)
# so the effect of barycentring on the Viterbi search can be compared later.
#
# Runs entirely on the host (PRESTO is a module, not in the container).
# Sources the project setup script for the environment.
#
# NOTE on dedispersion: the filterbank channels are already COHERENTLY
# dedispersed (within-channel) at the cluster DM 24.4 at recording time.
# prepdata performs the INCOHERENT dedispersion: it aligns the 128 channels
# in time at the chosen DM and sums them into a 1-D .dat time series.
# ----------------------------------------------------------------------------

# --- paths (edit if needed) -------------------------------------------------
SETUP="/fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh"

EXP_DIR="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/47Tuc_blind_search_v1"
STAGE1_DIR="${EXP_DIR}/stage1_clean"
CLEANED_FIL="${STAGE1_DIR}/47Tuc_22UL_1of2_L_RFIcleaned_01.fil"

# --- parameters -------------------------------------------------------------
DM="24.356441"
DM_TAG="DM24.356441"           # used in folder and file names
OUT_ROOT="47Tuc_${DM_TAG}"     # prepdata output rootname prefix

# --- environment ------------------------------------------------------------
echo "Start (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Host:        $(hostname)"
echo "Job ID:      ${SLURM_JOB_ID:-N/A}"
echo "Input:       ${CLEANED_FIL}"
echo "DM:          ${DM}"

# Source the environment WITHOUT 'set -e' (Lmod returns non-zero on harmless
# module warnings, which would otherwise abort the job).
set +e
source "${SETUP}"
set -e
set -uo pipefail

# Sanity: confirm prepdata is available
if ! command -v prepdata >/dev/null 2>&1; then
    echo "ERROR: prepdata not found on PATH after sourcing setup." >&2
    exit 1
fi
echo "prepdata: $(command -v prepdata)"

# Sanity: confirm input exists
if [ ! -f "${CLEANED_FIL}" ]; then
    echo "ERROR: cleaned filterbank not found: ${CLEANED_FIL}" >&2
    exit 1
fi

# --- output directories -----------------------------------------------------
BARY_DIR="${EXP_DIR}/stage2_dedisp/${DM_TAG}_bary"
TOPO_DIR="${EXP_DIR}/stage2_dedisp/${DM_TAG}_topo"
mkdir -p "${BARY_DIR}" "${TOPO_DIR}"

# ----------------------------------------------------------------------------
# Barycentred run (prepdata default = barycentred)
# ----------------------------------------------------------------------------
echo ""
echo "=== Barycentred dedispersion ==="
cd "${BARY_DIR}"
prepdata \
    -dm "${DM}" \
    -o "${OUT_ROOT}_bary" \
    "${CLEANED_FIL}" \
    2>&1 | tee "${BARY_DIR}/prepdata.log"

echo ""
echo "Barycentred outputs:"
ls -lh "${BARY_DIR}/${OUT_ROOT}_bary".dat "${BARY_DIR}/${OUT_ROOT}_bary".inf

# ----------------------------------------------------------------------------
# Topocentric run (prepdata -nobary)
# ----------------------------------------------------------------------------
echo ""
echo "=== Topocentric dedispersion (-nobary) ==="
cd "${TOPO_DIR}"
prepdata \
    -nobary \
    -dm "${DM}" \
    -o "${OUT_ROOT}_topo" \
    "${CLEANED_FIL}" \
    2>&1 | tee "${TOPO_DIR}/prepdata.log"

echo ""
echo "Topocentric outputs:"
ls -lh "${TOPO_DIR}/${OUT_ROOT}_topo".dat "${TOPO_DIR}/${OUT_ROOT}_topo".inf

# ----------------------------------------------------------------------------
echo ""
echo "=== Stage 2 complete ==="
echo "End (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Barycentred: ${BARY_DIR}/${OUT_ROOT}_bary.dat"
echo "Topocentric: ${TOPO_DIR}/${OUT_ROOT}_topo.dat"