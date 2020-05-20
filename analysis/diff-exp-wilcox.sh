#!/bin/bash
#SBATCH --job-name=diff-exp-wilcox
#SBATCH --output=diff-exp-wilcox.out
#SBATCH --error=diff-exp-wilcox.err
#SBATCH --time=24:00:00
#SBATCH --partition=gilad
#SBATCH --nodelist=midway-l16b-31
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=25600

module load R/3.4.3
R CMD BATCH --no-save --no-restore diff-exp-wilcox.R diff-exp-wilcox.out
