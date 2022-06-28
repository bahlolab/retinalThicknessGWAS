#!/bin/bash

#SBATCH -J adjustRefractiveError
#SBATCH -t 12:0:0
#SBATCH --mem=256G
#SBATCH --output=/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/logs/07_adjustRefractiveError.log


## load modules
module unload R
module load R/4.1.2

cd /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/scripts/

./07_adjustRefractiveError.R
