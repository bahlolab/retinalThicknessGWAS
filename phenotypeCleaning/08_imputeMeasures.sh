#!/bin/bash

#SBATCH -J imputePheno
#SBATCH -t 23:0:0
#SBATCH --mem=64G
#SBATCH --array 1-10
#SBATCH --output=/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/logs/08_imputeMeasures_%a.log


## load modules
module unload R
module load R/4.1.2

cd /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/scripts/

./08_imputeMeasures.R $SLURM_ARRAY_TASK_ID
