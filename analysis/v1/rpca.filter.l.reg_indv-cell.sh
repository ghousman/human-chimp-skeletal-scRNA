#!/bin/bash
#SBATCH --job-name=rpca.filter.l.reg_indv-cell
#SBATCH --output=rpca.filter.l.reg_indv-cell.out
#SBATCH --error=rpca.filter.l.reg_indv-cell.err
#SBATCH --time=03:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=6400

module load R/3.5.1
export PYTHONUSERBASE=$HOME/.local/python
R CMD BATCH --no-save --no-restore rpca.filter.l.reg_indv-cell.R rpca.filter.l.reg_indv-cell.out
