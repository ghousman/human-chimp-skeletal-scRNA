#!/bin/bash
#SBATCH --job-name=diff-exp-TEST
#SBATCH --output=diff-exp-TEST.out
#SBATCH --error=diff-exp-TEST.err
#SBATCH --time=15:00:00
#SBATCH --partition=bigmem2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=25600

module load R/3.4.3
R CMD BATCH --no-save --no-restore diff-exp-TEST.R diff-exp-TEST.out
