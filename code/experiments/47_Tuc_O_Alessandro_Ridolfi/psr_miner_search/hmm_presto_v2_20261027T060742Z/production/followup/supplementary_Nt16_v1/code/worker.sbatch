#!/bin/bash
set -o pipefail
SUPP_RUNNER=${SUPP_SRUN:-$(command -v srun)}
if [[ "${SUPP_MODE}" == native ]]; then
    source /fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh
    setup_status=$?
    if (( setup_status != 0 )); then exit "${setup_status}"; fi
    export PYTHONPATH="${PRESTO}/lib/python3.10/site-packages"
else
    module --force purge
    module load gcc/11.3.0
    module load python/3.10.4
fi
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 OMP_DYNAMIC=FALSE MPLBACKEND=Agg
export MPLCONFIGDIR="${SUPP_ROOT}/mpl_${SLURM_JOB_ID}"
SUPP_PYTHON=/fred/oz022/vdimarco/software/envrmnts/.venvs/viterbi_psr/bin/python
SUPP_INDEX=$(${SUPP_PYTHON} -c 'import json,os; print(json.load(open(os.environ["SUPP_BATCH"]))[int(os.environ.get("SLURM_ARRAY_TASK_ID",0))])')
"${SUPP_RUNNER}" --cpu-bind=cores "${SUPP_PYTHON}" "${SUPP_ROOT}/code/workflow.py" "${SUPP_MODE}" --index "${SUPP_INDEX}"
