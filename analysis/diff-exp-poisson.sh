#!/bin/bash
#SBATCH --job-name=diff-exp-poisson
#SBATCH --output=diff-exp-poisson.out
#SBATCH --error=diff-exp-poisson.err
#SBATCH --time=24:00:00
#SBATCH --partition=gilad
#SBATCH --nodelist=midway-l16b-31
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=25600

module load R/3.4.3
R CMD BATCH --no-save --no-restore diff-exp-poisson.R diff-exp-poisson.out
