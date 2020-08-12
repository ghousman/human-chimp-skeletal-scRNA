#!/bin/bash
#SBATCH --job-name=prep_spps-cell
#SBATCH --output=prep_spps-cell.out
#SBATCH --error=prep_spps-cell.err
#SBATCH --time=03:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=6400

module load R/3.5.1
export PYTHONUSERBASE=$HOME/.local/python
R CMD BATCH --no-save --no-restore prep_spps-cell.R prep_spps-cell.out
