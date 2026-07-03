#!/bin/bash
#SBATCH --job-name=T0_final_153
#SBATCH --output=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/47Tuc_blind_search_v1/inspection/cand_153_433hz/followup_v1/folds/job_outputs/%x_%A_%a.out
#SBATCH --error=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/47Tuc_blind_search_v1/inspection/cand_153_433hz/followup_v1/folds/job_outputs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=00:30:00
#SBATCH --array=0-13

set +e
source /fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh
set -e

T0_OFFSETS=(1500 1505 1510 1515 1520 1525 1530 1535 1540 1545 1550 1555 1560 1565)
T0_OFFSET=${T0_OFFSETS[$SLURM_ARRAY_TASK_ID]}

FIL="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/data/47Tuc_22UL_1of2_L.fil"
OUTDIR="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/47Tuc_blind_search_v1/inspection/cand_153_433hz/followup_v1/folds/T0_final"
mkdir -p ${OUTDIR}
cd ${OUTDIR}

T0=$(python3 -c "print(59391.2096286575 + ${T0_OFFSET}/86400.0)")

echo "Start (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "T0 offset: ${T0_OFFSET} s, T0 MJD: ${T0}"

prepfold \
    -f 433.358871357097 \
    -dm 25.33 \
    -bin \
    -pb 7117.938710 \
    -x 0.282 \
    -e 0.0 \
    -To ${T0} \
    -npart 64 -n 128 -nsub 64 \
    -nodmsearch -noxwin \
    -filterbank \
    -o cand153_T0final${T0_OFFSET} \
    ${FIL}

echo "End (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"