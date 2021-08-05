#!/bin/bash
sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=DEedgeR_${1}_${7}
#SBATCH --output=DEedgeR_${1}_${7}.out
#SBATCH --error=DEedgeR_${1}_${7}.err
#SBATCH --time=48:00:00
#SBATCH --partition=gilad
#SBATCH --nodelist=midway-l16b-31
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=83400

module load R/3.6.1
R CMD BATCH --no-save --no-restore '--args name="$1" filter.arg=$2 filter.type="$3" filter.param=$4 lfc=$5 fdr="$6" assign="$7"' DEedgeR.R DEedgeR_$1_$7.out

EOT


