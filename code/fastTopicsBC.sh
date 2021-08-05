#!/bin/bash
sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=fastTopicsBC_${1}_${2}
#SBATCH --output=fastTopicsBC_${1}_${2}.out
#SBATCH --error=fastTopicsBC_${1}_${2}.err
#SBATCH --time=96:00:00
#SBATCH --partition=gilad
#SBATCH --nodelist=midway-l16b-31
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=83400

module load R/3.6.1
R CMD BATCH --no-save --no-restore '--args clustnum=$1 data.source="$2"' fastTopicsBC.R fastTopicsBC_$1_$2.out

EOT


