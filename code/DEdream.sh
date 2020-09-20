#!/bin/bash
sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=DEdream_${1}_${8}
#SBATCH --output=DEdream_${1}_${8}.out
#SBATCH --error=DEdream_${1}_${8}.err
#SBATCH --time=48:00:00
#SBATCH --partition=gilad
#SBATCH --nodelist=midway-l16b-31
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=83400

module load R/3.6.1
R CMD BATCH --no-save --no-restore '--args name="$1" data.source="$2" data.type="$3" filter.arg=$4 filter.type="$5" filter.param=$6 fml="$7" lfc=FALSE fdr="BH" assign="$8"' DEdream.R DEdream_$1_$8.out

EOT


