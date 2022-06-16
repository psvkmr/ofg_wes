#!/bin/bash -l

#SBATCH --partition=brc,shared
#SBATCH --job-name=r_merge
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --output=/scratch/users/%u/tests/ofg_merge_chr_outputs_%J.out
#SBATCH --verbose

set -o noclobber

conda activate r4

Rscript /scratch/users/k2142172/scripts/misc_scripts/ofg_merge_chr_outputs.R
