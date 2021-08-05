#!/bin/bash
sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=FindMarkers_${1}
#SBATCH --output=FindMarkers_${1}.out
#SBATCH --error=FindMarkers_${1}.err
#SBATCH --time=36:00:00
#SBATCH --partition=gilad
#SBATCH --nodelist=midway-l16b-31
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=83400

module load R/3.6.1
R CMD BATCH --no-save --no-restore '--args res="$1"' FindMarkers.R FindMarkers_$1.out

EOT


