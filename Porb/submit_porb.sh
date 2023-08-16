#!/bin/bash
#
#SBATCH --job-name=porb
#SBATCH --output=logs/porb_%A_%a.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --tmp=10G
#SBATCH --array=0-49

source ~/activate_nt.sh
workon gc_search_nt

export PSRFITS_HEADER_DIR=/fred/oz002/ldunn/viterbi_psr_search/paper/Porb

mkdir -p $JOBFS/$SLURM_ARRAY_TASK_ID
cd $JOBFS/$SLURM_ARRAY_TASK_ID

cp /fred/oz002/ldunn/viterbi_psr_search/paper/Porb/noise_template.params .
cp /fred/oz002/ldunn/viterbi_psr_search/paper/Porb/sig_template.params .
cp /fred/oz002/ldunn/viterbi_psr_search/paper/Porb/template.par .

python /fred/oz002/ldunn/viterbi_psr_search/paper/Porb/scan_porb.py

cp porb_results.dat /fred/oz002/ldunn/viterbi_psr_search/paper/Porb/out/porb_results_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.dat
