#!/bin/bash
#SBATCH --job-name=de-edgeR
#SBATCH --output=de-edgeR.out
#SBATCH --error=de-edgeR.err
#SBATCH --time=36:00:00
#SBATCH --partition=gilad
#SBATCH --nodelist=midway-l16b-31
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=83400

module load R/3.4.3
R CMD BATCH --no-save --no-restore de-edgeR.R de-edgeR.out
