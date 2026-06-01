#!/bin/bash
#
#SBATCH --job-name=make_datasets
#SBATCH --output=logs/make_datasets_%A_%a.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --tmp=10G
#SBATCH --array=0-399

source ~/activate_nt.sh
workon gc_search_nt

export PSRFITS_HEADER_DIR=/fred/oz002/ldunn/viterbi_psr_search/paper/datasets

mkdir -p $JOBFS/$SLURM_ARRAY_TASK_ID
cd $JOBFS/$SLURM_ARRAY_TASK_ID

cp /fred/oz002/ldunn/viterbi_psr_search/paper/datasets/noise_template.params .
cp /fred/oz002/ldunn/viterbi_psr_search/paper/datasets/sig_template.params .
cp /fred/oz002/ldunn/viterbi_psr_search/paper/datasets/template.par .

python /fred/oz002/ldunn/viterbi_psr_search/paper/datasets/make_datasets.py $SLURM_ARRAY_TASK_ID

cp -r data_* /fred/oz002/ldunn/viterbi_psr_search/paper/datasets/datasets/
