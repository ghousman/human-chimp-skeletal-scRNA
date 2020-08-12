#!/bin/bash
#SBATCH --job-name=merge.l.6k
#SBATCH --output=merge.l.6k.out
#SBATCH --error=merge.l.6k.err
#SBATCH --time=03:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=6400

module load R/3.5.1
export PYTHONUSERBASE=$HOME/.local/python
R CMD BATCH --no-save --no-restore merge.l.6k.R merge.l.6k.out
