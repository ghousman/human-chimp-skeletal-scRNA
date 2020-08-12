#!/bin/bash
#SBATCH --job-name=diff-exp-mast
#SBATCH --output=diff-exp-mast.out
#SBATCH --error=diff-exp-mast.err
#SBATCH --time=24:00:00
#SBATCH --partition=gilad
#SBATCH --nodelist=midway-l16b-31
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=25600

module load R/3.6.1
R CMD BATCH --no-save --no-restore diff-exp-mast.R diff-exp-mast.out
