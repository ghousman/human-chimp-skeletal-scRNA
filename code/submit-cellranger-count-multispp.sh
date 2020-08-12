#!/bin/bash
sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=count_${1}_sbatch
#SBATCH --output=count_${1}_sbatch.out
#SBATCH --error=count_${1}_sbatch.err
#SBATCH --time=36:00:00
#SBATCH --partition=gilad
#SBATCH --nodelist=midway-l16b-31
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=20000

export PATH=/project2/gilad/ghousman/cellranger/cellranger-3.1.0:$PATH

cellranger count --id=$1 \
                 --transcriptome=$2 \
                 --fastqs=$3 \
		 --sample=$1 \
                 --expect-cells=5000
EOT
