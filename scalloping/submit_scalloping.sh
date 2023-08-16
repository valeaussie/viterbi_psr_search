#!/bin/bash
#
#SBATCH --job-name=scalloping
#SBATCH --output=logs/scalloping_%A_%a.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --tmp=10G
#SBATCH --array=0-49

source ~/activate_nt.sh
workon gc_search_nt

export PSRFITS_HEADER_DIR=/fred/oz002/ldunn/viterbi_psr_search/paper/scalloping

mkdir -p $JOBFS/$SLURM_ARRAY_TASK_ID
cd $JOBFS/$SLURM_ARRAY_TASK_ID

cp /fred/oz002/ldunn/viterbi_psr_search/paper/scalloping/noise_template.params .
cp /fred/oz002/ldunn/viterbi_psr_search/paper/scalloping/sig_template.params .
cp /fred/oz002/ldunn/viterbi_psr_search/paper/scalloping/template.par .

python /fred/oz002/ldunn/viterbi_psr_search/paper/scalloping/scan_scalloping.py

cp scalloping_results.dat /fred/oz002/ldunn/viterbi_psr_search/paper/scalloping/no_padding/out/scalloping_results_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.dat
