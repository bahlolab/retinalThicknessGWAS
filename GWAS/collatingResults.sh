#!/bin/bash

#SBATCH -J extractGWsig
#SBATCH -o /vast/scratch/users/jackson.v/retThickness/GWAS/logs/collatingResults_%A_%a.log
#SBATCH -t 10:0:0
#SBATCH --mem=100GB
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-22


chr=$SLURM_ARRAY_TASK_ID

module load R/4.1.3

./scripts/collatingResults.R  \
  --chr $chr 
