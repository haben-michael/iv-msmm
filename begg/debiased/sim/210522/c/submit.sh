#!/bin/bash
#SBATCH --time=3
#SBATCH --nodes=1
#SBATCH --mem=1GB
# use the current directory
#$ -cwd
#$ -S /bin/bash

#srun /home/users/habnice/R/bin/Rscript test.R
# srun module load R
# srun /home/users/habnice/.local/bin/Rscript $*
/usr/bin/Rscript $*


## submit command
# for b in {1..10}; do sbatch ./submit.sh sim.R 5e3; done; 
