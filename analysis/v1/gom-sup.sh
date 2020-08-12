#!/bin/bash
#SBATCH --job-name=gom-sup
#SBATCH --output=gom-sup.out
#SBATCH --error=gom-sup.err
#SBATCH --time=08:00:00
#SBATCH --partition=bigmem2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=51200

module load R/3.4.3
R CMD BATCH --no-save --no-restore gom-sup.R gom-sup.out
