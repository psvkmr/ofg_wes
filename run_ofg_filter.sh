#!/bin/bash -l

#SBATCH --partition=brc,shared
#SBATCH --job-name=r_filt
#SBATCH --time=04:00:00
#SBATCH --mem=96G
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=8
#SBATCH --output=/scratch/users/%u/tests/ofg_filter_chr_outputs_%J.out
#SBATCH --verbose

conda activate r4

Rscript /scratch/users/k2142172/scripts/misc_scripts/ofg_filter_chr_outputs.R
