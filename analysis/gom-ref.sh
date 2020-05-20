#!/bin/bash
#SBATCH --job-name=gom-ref
#SBATCH --output=gom-ref.out
#SBATCH --error=gom-ref.err
#SBATCH --time=12:00:00
#SBATCH --partition=bigmem2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=51200

module load R/3.4.3
R CMD BATCH --no-save --no-restore gom-ref.R gom-ref.out
