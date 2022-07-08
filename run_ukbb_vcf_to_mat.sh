#!/bin/bash

#SBATCH --partition=brc,shared
#SBATCH --job-name=vcfToMat
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --verbose
#SBATCH --output=/scratch/users/k2142172/tests/vcfToMat.out

Rscript ukbb_vcf_to_mat.R
