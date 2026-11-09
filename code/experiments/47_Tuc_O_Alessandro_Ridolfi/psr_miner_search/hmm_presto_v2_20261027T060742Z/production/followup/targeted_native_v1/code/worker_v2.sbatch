#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
set -o pipefail
TARGET_SRUN=$(command -v srun)
source /fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh
setup_status=$?
if (( setup_status != 0 )); then exit "${setup_status}"; fi
set -euo pipefail
export PYTHONPATH="${PRESTO}/lib/python3.10/site-packages"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg
args=()
if [[ -n "${TARGET_INDEX:-}" ]]; then args+=(--index "${TARGET_INDEX}"); fi
if [[ -n "${TARGET_BATCH:-}" ]]; then
  TARGET_INDEX=$(python -c 'import json,os; print(json.load(open(os.environ["TARGET_BATCH"]))[int(os.environ["SLURM_ARRAY_TASK_ID"])])')
  args+=(--index "${TARGET_INDEX}")
fi
if [[ -n "${TARGET_NAME:-}" ]]; then args+=(--target "${TARGET_NAME}"); fi
"${TARGET_SRUN}" --cpu-bind=cores python "${TARGET_CODE}/targeted_workflow_v2.py" "${TARGET_MODE}" "${args[@]}"
