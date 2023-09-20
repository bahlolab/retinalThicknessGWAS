#!/bin/bash

#SBATCH -J createSNPSummaryFile.sh
#SBATCH -o /vast/scratch/users/jackson.v/retThickness/GWAS/annot/error/createSummary_%A_%a.log
#SBATCH -t 8:0:0
#SBATCH --mem=36GB
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-23


chr=$SLURM_ARRAY_TASK_ID

module load R

cd /vast/scratch/users/jackson.v/retThickness/GWAS/annot/

./createSNPSummaryFile.R -c $chr