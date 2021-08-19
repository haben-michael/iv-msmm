#!/bin/bash
#SBATCH --time=3
#SBATCH --nodes=1
#SBATCH --mem=1GB
#$ -cwd
#$ -S /bin/bash

srun module load R
srun /home/users/$USER/.local/bin/Rscript $*

# example call on a shell with slurm installed
# for in {1..1}; do for n in $(seq 50 50 1000);do ./submit.sh sim.R 3 $n 1000 & done; done
