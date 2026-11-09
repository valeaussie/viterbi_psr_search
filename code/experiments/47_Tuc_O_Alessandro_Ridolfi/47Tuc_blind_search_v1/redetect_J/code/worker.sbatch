#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
J_SLURM_BIN=$(dirname "$(command -v srun)")
set +e
source /fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="${PRESTO}/lib/python3.10/site-packages:${J_CODE}:${PYTHONPATH:-}"
export PATH="${J_SLURM_BIN}:${PATH}"
export MPLCONFIGDIR="${JOBFS:-/tmp}/J-mpl-${SLURM_JOB_ID}"
export SLURM_TMPDIR="${JOBFS:-/tmp}"
J_INDEX=$(python -c 'import json,os;print(json.load(open(os.environ["J_BATCH"]))[int(os.environ["SLURM_ARRAY_TASK_ID"])])')
srun --cpu-bind=cores python "${J_CODE}/workflow.py" "${J_MODE}" --index "${J_INDEX}"
