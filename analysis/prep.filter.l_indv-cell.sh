#!/bin/bash
#SBATCH --job-name=prep.filter.l_indv-cell
#SBATCH --output=prep.filter.l_indv-cell.out
#SBATCH --error=prep.filter.l_indv-cell.err
#SBATCH --time=03:00:00
#SBATCH --partition=bigmem2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=25600

module load R/3.5.1
export PYTHONUSERBASE=$HOME/.local/python
R CMD BATCH --no-save --no-restore prep.filter.l_indv-cell.R prep.filter.l_indv-cell.out
