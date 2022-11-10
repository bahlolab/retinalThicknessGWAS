#!/bin/bash

#SBATCH -J fpcaAllCleanedScans
#SBATCH -t 47:59:0
#SBATCH --mem=256G
#SBATCH --output=/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenoExploratory/working/fpcaAllCleanedScans_var50Removed_%A.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jackson.v@wehi.edu.au


module load R/4.2.0

cd /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/retinalThicknessGWAS/generateFPCs/

./fpcaAllCleanedScans_var50Removed.R
