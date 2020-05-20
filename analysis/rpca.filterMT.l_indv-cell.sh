#!/bin/bash
#SBATCH --job-name=rpca.filterMT.l_indv-cell
#SBATCH --output=rpca.filterMT.l_indv-cell.out
#SBATCH --error=rpca.filterMT.l_indv-cell.err
#SBATCH --time=03:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=6400

module load R/3.5.1
export PYTHONUSERBASE=$HOME/.local/python
R CMD BATCH --no-save --no-restore rpca.filterMT.l_indv-cell.R rpca.filterMT.l_indv-cell.out
