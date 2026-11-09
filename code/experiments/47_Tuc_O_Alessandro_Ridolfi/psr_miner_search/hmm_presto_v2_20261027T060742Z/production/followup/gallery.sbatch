#!/bin/bash
#SBATCH --job-name=47tuc-track-gallery
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
set -o pipefail
set +e
source /fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh
setup_status=$?
set -euo pipefail
if (( setup_status != 0 )); then exit "${setup_status}"; fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg
python "${FOLLOWUP_ROOT}/recover_tracks.py" gallery --manifest "${FOLLOWUP_ROOT}/track_manifest.json"
