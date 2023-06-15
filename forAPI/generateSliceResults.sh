#!/bin/bash

#SBATCH -J formatSliceResults
#SBATCH -t 47:59:0
#SBATCH --mem=120G
#SBATCH --output=/vast/scratch/users/jackson.v/retThickness/GWAS/errors/formattingResultsYue_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-23


if [ "$SLURM_ARRAY_TASK_ID" = "23" ]
then
  CHR="X"
else
  CHR="$SLURM_ARRAY_TASK_ID"
fi

echo "CHR: $CHR"

module unload R
module load R/4.1.3

cd /vast/scratch/users/jackson.v/retThickness/GWAS

./generateSliceResults.R --chr $CHR

cd /vast/scratch/users/jackson.v/retThickness/GWAS/forYue/
tar -czvf chr$CHR.tar.gz chr$CHR/

rsync -av chr$CHR.tar.gz /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/finalResultsEUR

