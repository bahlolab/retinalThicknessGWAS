#!/bin/bash

#SBATCH -J imputePhenoUnadjusted
#SBATCH -t 6:0:0
#SBATCH --mem=24G
#SBATCH --array 1-46
#SBATCH --output=/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/logs/08_imputeMeasuresUnadjusted_%a.log


## load modules
module unload R
module load R/4.1.2

cd /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/scripts/

./08_imputeMeasuresUnadjusted.R $SLURM_ARRAY_TASK_ID
