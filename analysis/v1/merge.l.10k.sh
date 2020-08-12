#!/bin/bash
#SBATCH --job-name=merge.l.10k
#SBATCH --output=merge.l.10k.out
#SBATCH --error=merge.l.10k.err
#SBATCH --time=03:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=6400

module load R/3.5.1
export PYTHONUSERBASE=$HOME/.local/python
R CMD BATCH --no-save --no-restore merge.l.10k.R merge.l.10k.out
