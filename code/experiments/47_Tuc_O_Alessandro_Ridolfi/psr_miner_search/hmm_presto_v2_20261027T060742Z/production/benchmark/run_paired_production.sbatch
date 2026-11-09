#!/bin/bash
#SBATCH --job-name=47tuc-presto-v2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --chdir=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/psr_miner_search
#SBATCH --output=controlled-benchmark-%j.out

set -eo pipefail
BENCHMARK_DIR=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/psr_miner_search
EXPERIMENT_ROOT=${BENCHMARK_DIR%/psr_miner_search}
DATA_ROOT="$EXPERIMENT_ROOT/47Tuc_blind_search_v1/stage2_dedisp"
PACKAGE="${SEARCH_PACKAGE:?Set the frozen package directory}/benchmark"
RUN_ROOT="${SEARCH_PACKAGE}/paired_${SLURM_JOB_ID}"

date -u
hostname
BENCHMARK_SRUN=$(command -v srun)
set +e
source /fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh
setup_status=$?
set -e
if [ "$setup_status" -ne 0 ]; then exit "$setup_status"; fi
set -u
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLBACKEND=Agg
export MPLCONFIGDIR="$RUN_ROOT/matplotlib_cache"
"$BENCHMARK_SRUN" --cpu-bind=cores python "$PACKAGE/controlled_benchmark.py" hmm --root "$RUN_ROOT" --repeats 2

module load singularity
export APPTAINERENV_OMP_NUM_THREADS=1 APPTAINERENV_OMP_DYNAMIC=FALSE
export APPTAINERENV_OPENBLAS_NUM_THREADS=1 APPTAINERENV_MKL_NUM_THREADS=1 APPTAINERENV_NUMEXPR_NUM_THREADS=1
export APPTAINERENV_SLURM_JOB_ID="$SLURM_JOB_ID" APPTAINERENV_SLURM_CPUS_PER_TASK=1
SIF=/fred/oz005/users/ecarli/singularity_images/pulsar-miner-2.0-prestozl-22feb2025_latest.sif
stat "$SIF"
"$BENCHMARK_SRUN" --cpu-bind=cores singularity exec --cleanenv \
    -B "$EXPERIMENT_ROOT:$EXPERIMENT_ROOT" -B "$DATA_ROOT:$DATA_ROOT:ro" \
    "$SIF" python3 "$PACKAGE/controlled_benchmark.py" presto --root "$RUN_ROOT" --all-dms
date -u

python3 "$PACKAGE/verify_paired.py" "$RUN_ROOT" --expected-dms 21 --expected-repeats 2
