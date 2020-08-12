#!/bin/bash
#SBATCH --job-name=de-MAST
#SBATCH --output=de-MAST.out
#SBATCH --error=de-MAST.err
#SBATCH --time=24:00:00
#SBATCH --partition=gilad
#SBATCH --nodelist=midway-l16b-31
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=38400

module load R/3.6.1
R CMD BATCH --no-save --no-restore de-MAST.R de-MAST.out
