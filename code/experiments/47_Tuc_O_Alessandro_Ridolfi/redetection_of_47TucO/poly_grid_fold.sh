#!/bin/bash
#SBATCH --job-name=poly_grid_fold
#SBATCH --output=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/slurm_outputs/slurm-polyfold-%A_%a.out
#SBATCH --error=/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/slurm_outputs/slurm-polyfold-%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8g
#SBATCH --time=00:40:00
#SBATCH --array=0-14

set -o pipefail

CODE_DIR="/fred/oz022/vdimarco/software/install/viterbi_psr_search/code"
EXP="${CODE_DIR}/experiments/47_Tuc_O_Alessandro_Ridolfi"
REDET_DIR="${EXP}/redetection_of_47TucO"
PULSARX_SIF="${EXP}/data/pulsarx_20241219.sif"
BIND_PATH="/fred/oz022/vdimarco"
FIL="${EXP}/data/47Tuc_22UL_1of2_L.fil"
TEMPLATE="/home/pulsarx/software/PulsarX/include/template/meerkat_fold.template"
DM="24.356441"

ml singularity
cd "${REDET_DIR}"

NTS=(16 20 24 28 30 31 32 40 48 56 64 128 256 512 1024)
NT=${NTS[${SLURM_ARRAY_TASK_ID}]}

INFO="47TucO_redetect_Nt${NT}_psrfold.candfile.info"
[ -f "${INFO}" ] || { echo "missing ${INFO}"; exit 1; }

# Row 0 is the polynomial block; read its flags directly (no conversion).
PEPOCH=$(awk '/Row 0/{f=1} /Row 1/{f=0} f&&/--pepoch/{print $3}' "${INFO}")
F0=$(awk '/Row 0/{f=1} /Row 1/{f=0} f&&/--f0/{print $3}' "${INFO}")
F1=$(awk '/Row 0/{f=1} /Row 1/{f=0} f&&/--f1/{print $3}' "${INFO}")
F2=$(awk '/Row 0/{f=1} /Row 1/{f=0} f&&/--f2/{print $3}' "${INFO}")
echo "Nt=${NT} pepoch=${PEPOCH} f0=${F0} f1=${F1} f2=${F2}"

singularity exec -B "${BIND_PATH}" "${PULSARX_SIF}" psrfold_fil \
  -v -t 2 \
  --f0 "${F0}" --f1 "${F1}" --f2 "${F2}" \
  --dm "${DM}" --pepoch "${PEPOCH}" \
  --clfd 2.0 --rfi kadaneF 8 4 zdot --fill mean --fillPatch mean \
  -n 256 -b 64 -L 30 \
  --template "${TEMPLATE}" \
  --srcname polynomial_Nt${NT} \
  -o 47TucO_polynomial_Nt${NT} \
  -f "${FIL}" \
  2>&1 | tee 47TucO_polynomial_Nt${NT}.log