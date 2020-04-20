#!/bin/bash
#SBATCH --time=10
#SBATCH --nodes=1
#SBATCH --mem=4GB
#$ -cwd
#$ -S /bin/bash

\srun ~/.local/bin/Rscript $*


## submit command
# for n in $(seq 30 50 230); do for b in {1..10}; do sbatch ./submit.sh sim.R $n; done; done
