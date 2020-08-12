#!/bin/bash
#SBATCH --job-name=diff-exp-wilcox0.1
#SBATCH --output=diff-exp-wilcox0.1.out
#SBATCH --error=diff-exp-wilcox0.1.err
#SBATCH --time=24:00:00
#SBATCH --partition=gilad
#SBATCH --nodelist=midway-l16b-31
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=25600

module load R/3.4.3
R CMD BATCH --no-save --no-restore diff-exp-wilcox.R diff-exp-wilcox0.1.out
