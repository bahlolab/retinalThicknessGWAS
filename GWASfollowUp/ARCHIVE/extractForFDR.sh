#!/bin/bash

#SBATCH -J extractForFDR
#SBATCH -o /vast/scratch/users/jackson.v/retThickness/GWAS/logs/extractForFDR_%A_%a.log
#SBATCH -t 6:0:0
#SBATCH --mem=500MB
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-23

if [[ $SLURM_ARRAY_TASK_ID -eq 23 ]]; then
  CHR=X
else
  CHR=$SLURM_ARRAY_TASK_ID
fi

cd /vast/scratch/users/jackson.v/retThickness/GWAS/FDR

## manual test
# CHR=1
# SLICE=1
# PIX=1_50
# file=/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/results/chr${CHR}/${SLICE}/chr${CHR}Pixel.${PIX}.glm.linear.gz


# zcat  $file | awk -v var=$PIX '$7 < 5E-5 {printf "%s %s %s\n", $2, $7, var}'  > chr${CHR}test.txt

while read PIX SLICE y ; do
  echo "pix=$PIX"

file=/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/results/chr${CHR}/${SLICE}/chr${CHR}Pixel.${PIX}.glm.linear.gz
zcat  $file | awk -v var=$PIX '$7 < 5E-5 {printf "%s %s %s\n", $2, $7, var}'  >> chr${CHR}.txt

done < pixels.txt
