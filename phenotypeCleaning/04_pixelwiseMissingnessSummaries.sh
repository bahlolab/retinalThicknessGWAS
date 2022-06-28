#!/bin/bash

#SBATCH -J pixelwiseMissingnessSummaries
#SBATCH -t 2:0:0
#SBATCH --mem=96G
#SBATCH --output=/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/logs/04_pixelwiseMissingnessSummaries.log


## load modules
module unload R
module load R/4.1.2

cd /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/scripts/

./04_pixelwiseMissingnessSummaries.R
