#!/bin/bash
#SBATCH --job-name=47tuc-track-recovery
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
set -o pipefail
RECOVERY_SRUN=$(command -v srun)
set +e
source /fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh
setup_status=$?
set -euo pipefail
if (( setup_status != 0 )); then exit "${setup_status}"; fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 OMP_DYNAMIC=FALSE MPLBACKEND=Agg
export MPLCONFIGDIR="${FOLLOWUP_ROOT}/mpl_cache_${SLURM_ARRAY_TASK_ID}"
"${RECOVERY_SRUN}" --cpu-bind=cores python "${FOLLOWUP_ROOT}/recover_tracks.py" recover --manifest "${FOLLOWUP_ROOT}/track_manifest.json" --index "${SLURM_ARRAY_TASK_ID}"
