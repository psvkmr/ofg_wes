#!/bin/bash -l

#SBATCH --partition=brc,shared
#SBATCH --job-name=vcfToMat
#SBATCH --time=04:00:00
#SBATCH --mem=48G
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --verbose
#SBATCH --output=/scratch/users/k2142172/tests/vcfToMat.out

conda activate r4

Rscript /scratch/users/k2142172/scripts/ofg_wes/ukbb_vcf_to_mat.R
