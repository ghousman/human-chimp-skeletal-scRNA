#!/bin/bash
#SBATCH --job-name=rpca_pair-cell
#SBATCH --output=rpca_pair-cell.out
#SBATCH --error=rpca_pair-cell.err
#SBATCH --time=03:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=6400

module load R/3.5.1
export PYTHONUSERBASE=$HOME/.local/python
R CMD BATCH --no-save --no-restore rpca_pair-cell.R rpca_pair-cell.out
