#!/bin/bash
#SBATCH --job-name=corrmatrix
#SBATCH --output=corrmatrix.out
#SBATCH --error=corrmatrix.err
#SBATCH --time=04:00:00
#SBATCH --partition=bigmem2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=12800

module load R/3.4.3
R CMD BATCH --no-save --no-restore corrmatrix.R corrmatrix.out
