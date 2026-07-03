#!/bin/bash
#SBATCH --job-name=T0_grid_163
#SBATCH --output=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/47Tuc_blind_search_v1/inspection/cand_163_650hz/followup_v1/folds/job_outputs/%x_%A_%a.out
#SBATCH --error=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/47Tuc_blind_search_v1/inspection/cand_163_650hz/followup_v1/folds/job_outputs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=00:30:00
#SBATCH --array=0-20

set +e
source /fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh
set -e

# Grid over one full orbit (21420s) in steps of 1000s
T0_OFFSETS=(0 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 11000 12000 13000 14000 15000 16000 17000 18000 19000 20000)
T0_OFFSET=${T0_OFFSETS[$SLURM_ARRAY_TASK_ID]}

FIL="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/data/47Tuc_22UL_1of2_L.fil"
OUTDIR="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/47Tuc_blind_search_v1/inspection/cand_163_650hz/followup_v1/folds/T0_grid"
mkdir -p ${OUTDIR}
cd ${OUTDIR}

T0=$(python3 -c "print(59391.2156 + ${T0_OFFSET}/86400.0)")

echo "Start (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "T0 offset: ${T0_OFFSET} s, T0 MJD: ${T0}"

prepfold \
    -f 649.754437119645 \
    -dm 25.20 \
    -bin \
    -pb 21419.816208 \
    -x 1.33 \
    -e 0.0 \
    -To ${T0} \
    -npart 64 -n 128 -nsub 64 \
    -nodmsearch -noxwin \
    -filterbank \
    -o cand163_T0offset${T0_OFFSET} \
    ${FIL}

echo "End (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"