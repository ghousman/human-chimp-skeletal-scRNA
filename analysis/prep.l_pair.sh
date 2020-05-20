#!/bin/bash
#SBATCH --job-name=prep.l_pair
#SBATCH --output=prep.l_pair.out
#SBATCH --error=prep.l_pair.err
#SBATCH --time=03:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=6400

module load R/3.5.1
export PYTHONUSERBASE=$HOME/.local/python
R CMD BATCH --no-save --no-restore prep.l_pair.R prep.l_pair.out
