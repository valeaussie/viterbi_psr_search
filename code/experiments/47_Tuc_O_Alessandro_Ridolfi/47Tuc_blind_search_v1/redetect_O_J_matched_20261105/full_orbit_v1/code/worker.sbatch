#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --tmp=1G
O_FULL_SLURM_BIN=$(dirname "$(command -v srun)")
source /fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh
O_FULL_SETUP_STATUS=$?
if (( O_FULL_SETUP_STATUS != 0 )); then exit "${O_FULL_SETUP_STATUS}"; fi
set -euo pipefail
export PATH="${O_FULL_SLURM_BIN}:${PATH}"
export PYTHONPATH="${PRESTO}/lib/python3.10/site-packages:${J_FULL_CODE}:${PYTHONPATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1
export MPLCONFIGDIR="${JOBFS}/O_full_mpl"
if [[ "${J_FULL_MODE}" == "build" ]]; then
    srun --cpu-bind=cores python "${J_FULL_CODE}/full_orbit.py" prepare
    srun --cpu-bind=cores python "${J_FULL_CODE}/full_orbit.py" build
elif [[ "${J_FULL_MODE}" == "spin" ]]; then
    srun --cpu-bind=cores python "${J_FULL_CODE}/spin_prepare.py"
    srun --cpu-bind=cores python "${J_FULL_CODE}/../spin_refined/code/full_orbit.py" build
    srun --cpu-bind=cores python "${J_FULL_CODE}/../spin_refined/code/full_orbit.py" frequency_refinement
elif [[ "${J_FULL_MODE}" == "report" ]]; then
    srun --cpu-bind=cores python "${J_FULL_CODE}/report_all.py"
else
    srun --cpu-bind=cores python "${J_FULL_CODE}/full_orbit.py" "${J_FULL_MODE}"
fi
