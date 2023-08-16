#!/bin/bash
#
#SBATCH --job-name=flux
#SBATCH --output=logs/flux_%A_%a.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --tmp=10G
#SBATCH --array=0-49

source ~/activate_nt.sh
workon gc_search_nt

export PSRFITS_HEADER_DIR=/fred/oz002/ldunn/viterbi_psr_search/paper/flux

mkdir -p $JOBFS/$SLURM_ARRAY_TASK_ID
cd $JOBFS/$SLURM_ARRAY_TASK_ID

cp /fred/oz002/ldunn/viterbi_psr_search/paper/flux/noise_template.params .
cp /fred/oz002/ldunn/viterbi_psr_search/paper/flux/sig_template.params .
cp /fred/oz002/ldunn/viterbi_psr_search/paper/flux/template.par .

python /fred/oz002/ldunn/viterbi_psr_search/paper/flux/scan_flux.py

cp flux_results.dat /fred/oz002/ldunn/viterbi_psr_search/paper/flux/out/flux_results_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.dat
