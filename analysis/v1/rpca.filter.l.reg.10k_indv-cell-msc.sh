#!/bin/bash
#SBATCH --job-name=rpca.filter.l.reg.10k_indv-cell-msc
#SBATCH --output=rpca.filter.l.reg.10k_indv-cell-msc.out
#SBATCH --error=rpca.filter.l.reg.10k_indv-cell-msc.err
#SBATCH --time=12:00:00
#SBATCH --partition=bigmem2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=12800

module load R/3.5.1
export PYTHONUSERBASE=$HOME/.local/python
R CMD BATCH --no-save --no-restore rpca.filter.l.reg.10k_indv-cell-msc.R rpca.filter.l.reg.10k_indv-cell-msc.out
