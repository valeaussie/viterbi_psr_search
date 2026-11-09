#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
set -o pipefail
CANDIDATE_SRUN=$(command -v srun)
source /fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh
setup_status=$?
if (( setup_status != 0 )); then exit "${setup_status}"; fi
set -euo pipefail
export PYTHONPATH="${PRESTO}/lib/python3.10/site-packages"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg
if [[ -n "${CANDIDATE_BATCH:-}" ]]; then
  CANDIDATE_INDEX=$(python -c 'import json,os;print(json.load(open(os.environ["CANDIDATE_BATCH"]))[int(os.environ["SLURM_ARRAY_TASK_ID"])])')
fi
candidate_args=()
if [[ "${CANDIDATE_PILOT:-0}" == 1 ]]; then candidate_args+=(--pilot); fi
"${CANDIDATE_SRUN}" --cpu-bind=cores python "${CANDIDATE_CODE}/candidate_folds_v3.py" fold --index "${CANDIDATE_INDEX}" "${candidate_args[@]}"
