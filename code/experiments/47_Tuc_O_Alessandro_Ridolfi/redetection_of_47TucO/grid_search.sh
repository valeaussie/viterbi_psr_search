#!/bin/bash
#SBATCH --job-name=grid_search
#SBATCH --output=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/slurm_outputs/slurm-gridsearch-%A_%a.out
#SBATCH --error=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/slurm_outputs/slurm-gridsearch-%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8g
#SBATCH --time=00:40:00
#SBATCH --array=0-14

set -o pipefail

SETUP="/fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh"
CODE_DIR="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code"
EXP="${CODE_DIR}/experiments/47_Tuc_O_Alessandro_Ridolfi"
REDET_DIR="${EXP}/redetection_of_47TucO"
SEARCH="${CODE_DIR}/viterbi_search_and_fit.py"
DAT="47TucO_redetect_cfbf00000_Plan1_1_DM24.36.dat"

# Fixed observation start epoch (from the DM24.36 .inf).
MJD0="59391.125996527684038"

source "${SETUP}"
cd "${REDET_DIR}"

NTS=(16 20 24 28 30 31 32 40 48 56 64 128 256 512 1024)
TSFTS=(899.713885 719.771093 599.809282 514.122187 479.847395 464.368420 449.856943 359.885508 299.904603 257.061055 224.928433 112.464217 56.232108 28.116016 14.058008)

NT=${NTS[${SLURM_ARRAY_TASK_ID}]}
TSFT=${TSFTS[${SLURM_ARRAY_TASK_ID}]}
echo "Nt=${NT} Tsft=${TSFT} mjd0=${MJD0}"

python3 "${SEARCH}" \
  --infile "${DAT}" \
  --tsamp 0.00007656074766355 \
  --Tsft ${TSFT} --Nsft ${NT} \
  --f0 376.31 --bw 4.0 --poly-order 5 --plot-path \
  --dm 24.356441 --mjd-start ${MJD0} \
  --known-label "J0024-7204O" \
  --out_prefix 47TucO_redetect_Nt${NT} \
  2>&1 | tee 47TucO_redetect_Nt${NT}.log