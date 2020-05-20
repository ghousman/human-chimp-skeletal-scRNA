#!/bin/bash
#SBATCH --job-name=prep.filter.l_indv
#SBATCH --output=prep.filter.l_indv.out
#SBATCH --error=prep.filter.l_indv.err
#SBATCH --time=03:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=6400

module load R/3.5.1
export PYTHONUSERBASE=$HOME/.local/python
R CMD BATCH --no-save --no-restore prep.filter.l_indv.R prep.filter.l_indv.out
