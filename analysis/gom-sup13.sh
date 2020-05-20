#!/bin/bash
#SBATCH --job-name=gom-sup13
#SBATCH --output=gom-sup13.out
#SBATCH --error=gom-sup13.err
#SBATCH --time=12:00:00
#SBATCH --partition=bigmem2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=51200

module load R/3.4.3
R CMD BATCH --no-save --no-restore gom-sup13.R gom-sup13.out
