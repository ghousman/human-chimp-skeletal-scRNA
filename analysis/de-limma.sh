#!/bin/bash
#SBATCH --job-name=de-limma
#SBATCH --output=de-limma.out
#SBATCH --error=de-limma.err
#SBATCH --time=24:00:00
#SBATCH --partition=gilad
#SBATCH --nodelist=midway-l16b-31
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=38400

module load R/3.4.3
R CMD BATCH --no-save --no-restore de-limma.R de-limma.out
