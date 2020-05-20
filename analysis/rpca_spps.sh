#!/bin/bash
#SBATCH --job-name=rpca_spps
#SBATCH --output=rpca_spps.out
#SBATCH --error=rpca_spps.err
#SBATCH --time=03:00:00
#SBATCH --partition=bigmem2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=12800

module load R/3.5.1
export PYTHONUSERBASE=$HOME/.local/python
R CMD BATCH --no-save --no-restore rpca_spps.R rpca_spps.out
