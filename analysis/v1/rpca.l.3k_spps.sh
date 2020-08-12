#!/bin/bash
#SBATCH --job-name=rpca.l.3k_spps
#SBATCH --output=rpca.l.3k_spps.out
#SBATCH --error=rpca.l.3k_spps.err
#SBATCH --time=03:00:00
#SBATCH --partition=bigmem2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=12800

module load R/3.5.1
export PYTHONUSERBASE=$HOME/.local/python
R CMD BATCH --no-save --no-restore rpca.l.3k_spps.R rpca.l.3k_spps.out
