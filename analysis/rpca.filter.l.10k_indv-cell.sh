#!/bin/bash
#SBATCH --job-name=rpca.filter.l.10k_indv-cell
#SBATCH --output=rpca.filter.l.10k_indv-cell.out
#SBATCH --error=rpca.filter.l.10k_indv-cell.err
#SBATCH --time=03:00:00
#SBATCH --partition=bigmem2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=12800

module load R/3.5.1
export PYTHONUSERBASE=$HOME/.local/python
R CMD BATCH --no-save --no-restore rpca.filter.l.10k_indv-cell.R rpca.filter.l.10k_indv-cell.out
