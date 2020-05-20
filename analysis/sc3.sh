#!/bin/bash
#SBATCH --job-name=sc3
#SBATCH --output=sc3.out
#SBATCH --error=sc3.err
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=6400

module load R/3.4.3
R CMD BATCH --no-save --no-restore sc3.R sc3.out
