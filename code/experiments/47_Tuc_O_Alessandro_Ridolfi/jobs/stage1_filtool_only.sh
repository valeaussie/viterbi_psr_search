#!/bin/bash
#SBATCH --job-name=stage1_filtool
#SBATCH --output=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/slurm_outputs/slurm-filtool-%j.out
#SBATCH --error=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/slurm_outputs/slurm-filtool-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16g
#SBATCH --tmp=40g
#SBATCH --time=02:00:00

# ----------------------------------------------------------------------------
# Stage 1 (filtool only): RFI cleaning with PulsarX filtool inside Apptainer.
# Nothing else. No Python, no post-processing.
# ----------------------------------------------------------------------------

set -euo pipefail

# --- paths (edit if needed) -------------------------------------------------
INPUT_FIL="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/data/47Tuc_22UL_1of2_L.fil"
OUT_DIR="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/47Tuc_blind_search_v1/stage1_clean"
PULSARX_SIF="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/data/pulsarx_20241219.sif"
BIND_PATH="/fred/oz022/vdimarco"
OUTPUT_ROOT="47Tuc_22UL_1of2_L_RFIcleaned"
NTHREADS=8

# RA / Dec hard-coded from the header (no Python needed)
RA="00:24:04.6500"
DEC="-72:04:53.8000"

RFI_FLAGS="-z zdot kadaneF 8 4 kadaneT 8 4 zap 925 960 zap 1525 1612 zap 1675 1720 --fillPatch rand"

# --- run --------------------------------------------------------------------
mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

echo "Start (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Host:        $(hostname)"
echo "Job ID:      ${SLURM_JOB_ID:-N/A}"

module load apptainer

apptainer exec -B "${BIND_PATH}" "${PULSARX_SIF}" \
    filtool \
        --threads "${NTHREADS}" \
        --ra "${RA}" \
        --dec "${DEC}" \
        --telescope MeerKAT \
        --rootname "${OUTPUT_ROOT}" \
        --source_name J0024-7204O \
        ${RFI_FLAGS} \
        -f "${INPUT_FIL}" \
        2>&1 | tee "${OUT_DIR}/filtool.log"
        
echo "" >> "${OUT_DIR}/filtool.log"
echo "# NOTE: 'nsamples : 26122' above is filtool's PER-SEGMENT count (seglen=2s)," >> "${OUT_DIR}/filtool.log"
echo "# not the total. Total nsamples in output ~ 188e6 (full 4-hour observation)." >> "${OUT_DIR}/filtool.log"
echo ""
echo "End (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"
ls -lh "${OUT_DIR}/${OUTPUT_ROOT}"_*.fil