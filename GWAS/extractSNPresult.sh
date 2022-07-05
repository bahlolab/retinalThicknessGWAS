#!/bin/bash

#SBATCH -J extractSNPResult
#SBATCH -o /vast/scratch/users/jackson.v/retThickness/GWAS/logs/extractResultrs17421627_%A_%a.log
#SBATCH -t 6:0:0
#SBATCH --mem=10MB
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-119


slice=$SLURM_ARRAY_TASK_ID
outDir=/vast/scratch/users/jackson.v/retThickness/GWAS/results/
chr=5
snp=rs17421627

resultsDir=/vast/scratch/users/jackson.v/retThickness/GWAS/results/chr${chr}

# pixel=64_127

echo " outputting results for $snp"

cd $outDir

mkdir -p ${snp}Result

rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/pixels.txt .

nPix=$(wc -l pixels.txt | cut -d " " -f 1)

echo -e "$(head -n1 ${resultsDir}/${slice}/chr${chr}Pixel.${slice}_64.glm.linear) \t pixel" > ${snp}Result/${snp}Result_${slice}.txt

for pix in $(seq 1 $nPix)
do

  read pixel y x < <(sed -n ${pix}p pixels.txt)

  if [ $y -eq $slice ]
  then
    echo 'Pixel '$pixel''

    file=${resultsDir}/${slice}/chr${chr}Pixel.${pixel}.glm.linear

   awk -v id="$pixel" -v snp="$snp" '$3 == snp { print $0 ,"\t", id }' $file >> ${snp}Result/${snp}Result_${slice}.txt
 fi

done
