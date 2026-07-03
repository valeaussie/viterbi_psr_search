#!/bin/bash
#SBATCH --job-name=viterbi_DMscan
#SBATCH --output=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/47Tuc_blind_search_v1/inspection/cand_153_433hz/followup_v1/DM0_check/job_outputs/%x_%A_%a.out
#SBATCH --error=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/47Tuc_blind_search_v1/inspection/cand_153_433hz/followup_v1/DM0_check/job_outputs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=00:30:00
#SBATCH --array=0-4

set +e
source /fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh
set -e

DMS=(5.00 10.00 15.00 20.00 25.00)
DM=${DMS[$SLURM_ARRAY_TASK_ID]}
DMTAG=$(printf "%.2f" ${DM})

CODE_DIR="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code"
EXP="${CODE_DIR}/experiments/47_Tuc_O_Alessandro_Ridolfi/47Tuc_blind_search_v1"
VITERBI="${CODE_DIR}/viterbi_pipeline.py"
DMCHECK="${EXP}/inspection/cand_153_433hz/followup_v1/DM0_check"

DAT="${DMCHECK}/J0000-00_cfbf00000_Plan1_1_DM${DMTAG}.dat"
OUT="${DMCHECK}/viterbi_DM${DMTAG}"

cat > ${DMCHECK}/viterbi_DM${DMTAG}.params << EOF
infile         = "${DAT}"
tsamp          = 7.656074766355e-05
Tsft           = 56.23203125
Nsft           = -1
f0             = 430.32
bw             = 6.0
padding_factor = 1
num_harm       = 1
top_paths      = 1
save_delta     = False
out_prefix     = "${OUT}"
plot_path      = True
dm             = ${DM}
mjd_start      = 59391.125996527684038
EOF

echo "Running Viterbi at DM=${DM}..."
cd ${DMCHECK}
python3 ${VITERBI} --params ${DMCHECK}/viterbi_DM${DMTAG}.params
echo "Done DM=${DM}"