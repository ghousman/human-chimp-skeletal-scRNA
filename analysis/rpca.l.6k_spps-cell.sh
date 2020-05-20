#!/bin/bash
#SBATCH --job-name=rpca.l.6k_spps-cell
#SBATCH --output=rpca.l.6k_spps-cell.out
#SBATCH --error=rpca.l.6k_spps-cell.err
#SBATCH --time=03:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=6400

module load R/3.5.1
export PYTHONUSERBASE=$HOME/.local/python
R CMD BATCH --no-save --no-restore rpca.l.6k_spps-cell.R rpca.l.6k_spps-cell.out
