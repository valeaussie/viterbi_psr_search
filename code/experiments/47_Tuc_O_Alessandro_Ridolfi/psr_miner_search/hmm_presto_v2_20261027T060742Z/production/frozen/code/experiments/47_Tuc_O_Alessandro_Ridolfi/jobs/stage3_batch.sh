#!/bin/bash
#SBATCH --job-name=47tuc-hmm-v2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
set -euo pipefail
: "${SEARCH_CODE_DIR:?Set the frozen source code directory}"
: "${SEARCH_OUTPUT_BASE:?Set the new output directory}"
: "${SEARCH_GRID_SIZE:?Set from the frozen grid manifest}"
N_PER_WORKER=${SEARCH_CONFIGS_PER_WORKER:-16}
FIRST=$(( SLURM_ARRAY_TASK_ID * N_PER_WORKER ))
LAST=$(( FIRST + N_PER_WORKER ))
if (( LAST > SEARCH_GRID_SIZE )); then LAST=${SEARCH_GRID_SIZE}; fi
SCRIPT="${SEARCH_CODE_DIR}/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/stage3_blind.sh"
for (( SEARCH_GRID_INDEX=FIRST; SEARCH_GRID_INDEX<LAST; SEARCH_GRID_INDEX++ )); do
    export SEARCH_GRID_INDEX
    bash "${SCRIPT}"
done
