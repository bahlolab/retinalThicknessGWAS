#!/bin/bash

#SBATCH -J extractGWsig
#SBATCH -o /vast/scratch/users/jackson.v/retThickness/GWAS/logs/extractResultGWsig_%A_%a.log
#SBATCH -t 1:0:0
#SBATCH --mem=10MB
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-119


slice=$SLURM_ARRAY_TASK_ID
outDir=/vast/scratch/users/jackson.v/retThickness/GWAS/GWsigResults/

for chr in {1..22}
do

# chr=5

resultsDir=/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/results/chr${chr}/$slice

echo " outputting results for chr$chr "

cd $outDir

mkdir -p chr${chr}

awk -v y="$slice" ' $2==y ' /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/pixels.txt > ${slice}_pixels.txt

nPix=$(wc -l  ${slice}_pixels.txt | cut -d " " -f 1)

echo -e "$(zcat ${resultsDir}/chr${chr}Pixel.${slice}_64.glm.linear.gz | head -n1) \t pixel" > chr${chr}/chr${chr}Slice${slice}_5e-5Sig.txt


for pix in $(seq 1 $nPix)
do

  read pixel y x < <(sed -n ${pix}p ${slice}_pixels.txt)

    file=${resultsDir}/chr${chr}Pixel.${pixel}.glm.linear.gz
    ref=/vast/scratch/users/jackson.v/retThickness/GWAS/GWsigResults/chr${chr}/chr${chr}all5e-5SigSNPs.txt
 
   awk -v id="$pixel" 'FNR==NR {f1[$1]; next} $2 in f1 { print $0 ,"\t", id }' $ref <(zcat $file)  >> chr${chr}/chr${chr}Slice${slice}_5e-5Sig.txt

done

done

chr=X

resultsDir=/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/results/chr${chr}/$slice

echo " outputting results for chr$chr "

cd $outDir

mkdir -p chr${chr}

awk -v y="$slice" ' $2==y ' /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/pixels.txt > ${slice}_pixels.txt

nPix=$(wc -l  ${slice}_pixels.txt | cut -d " " -f 1)

echo -e "$(zcat ${resultsDir}/chr${chr}Pixel.${slice}_64.glm.linear.gz | head -n1) \t pixel" > chr${chr}/chr${chr}Slice${slice}_5e-5Sig.txt


for pix in $(seq 1 $nPix)
do

  read pixel y x < <(sed -n ${pix}p ${slice}_pixels.txt)

    file=${resultsDir}/chr${chr}Pixel.${pixel}.glm.linear.gz
    ref=/vast/scratch/users/jackson.v/retThickness/GWAS/GWsigResults/chr${chr}/chr${chr}all5e-5SigSNPs.txt
 
   awk -v id="$pixel" 'FNR==NR {f1[$1]; next} $2 in f1 { print $0 ,"\t", id }' $ref <(zcat $file)  >> chr${chr}/chr${chr}Slice${slice}_5e-5Sig.txt

done
