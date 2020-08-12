#!/bin/bash
#SBATCH --job-name=corrmatrix3
#SBATCH --output=corrmatrix3.out
#SBATCH --error=corrmatrix3.err
#SBATCH --time=12:00:00
#SBATCH --partition=bigmem2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=12800

module load R/3.4.3
R CMD BATCH --no-save --no-restore corrmatrix3.R corrmatrix3.out
