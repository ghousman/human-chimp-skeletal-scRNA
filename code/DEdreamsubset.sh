#!/bin/bash
sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=DEdream_${1}_${9}
#SBATCH --output=DEdream_${1}_${9}.out
#SBATCH --error=DEdream_${1}_${9}.err
#SBATCH --time=48:00:00
#SBATCH --partition=gilad
#SBATCH --nodelist=midway-l16b-31
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=83400

module load R/3.6.1
R CMD BATCH --no-save --no-restore '--args name="$1" data.source="$2" subset.rep="$3" subset.num=$4 filter.arg=$5 filter.type="$6" filter.param=$7 fml="$8" lfc=FALSE fdr="BH" assign="$9"' DEdreamsubset.R DEdream_$1_$9.out

EOT


