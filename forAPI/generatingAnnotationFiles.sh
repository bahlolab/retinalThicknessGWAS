#!/bin/bash

workDir=/vast/scratch/users/jackson.v/retThickness/GWAS
dataDir=/vast/scratch/users/jackson.v/retThickness/GWAS/data/cleanedEURData

cd $workDir

 mkdir -p $workDir/plink
 mkdir -p $workDir/data
 cd $workDir/plink
 wget "https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_x86_64_20221024.zip"
 unzip plink2_linux_x86_64_20221024.zip

 mkdir -p /vast/scratch/users/jackson.v/retThickness/GWAS/annot
 mkdir -p $workDir/annot/data/

cd /vast/scratch/users/jackson.v/retThickness/GWAS/annot


for chr in {1..22} X
do

 ## all SNPs
../plink/plink2 \
  --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr${chr} \
  --make-just-pvar cols=vcfheader \
  --out chr${chr}

../plink/plink2 \
 --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr${chr} \
 --freq \
 --out chr${chr}Freq

done



