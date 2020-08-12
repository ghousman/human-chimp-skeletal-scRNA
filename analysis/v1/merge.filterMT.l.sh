#!/bin/bash
#SBATCH --job-name=merge.filterMT.l
#SBATCH --output=merge.filterMT.l.out
#SBATCH --error=merge.filterMT.l.err
#SBATCH --time=03:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=6400

module load R/3.5.1
export PYTHONUSERBASE=$HOME/.local/python
R CMD BATCH --no-save --no-restore merge.filterMT.l.R merge.filterMT.l.out
