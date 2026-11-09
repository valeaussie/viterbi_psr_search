#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
J_SLURM_BIN=$(dirname "$(command -v srun)")
set +e
source /fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh
set -euo pipefail
export PATH="${J_SLURM_BIN}:${PATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="${PRESTO}/lib/python3.10/site-packages:${PYTHONPATH:-}"
export MPLCONFIGDIR="${JOBFS:-/tmp}/J-full-mpl-${SLURM_JOB_ID}"
srun --cpu-bind=cores python "${J_FULL_CODE}/full_orbit.py" "${J_FULL_MODE}"
