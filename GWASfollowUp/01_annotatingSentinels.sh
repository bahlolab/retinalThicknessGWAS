#!/bin/bash

workDir=/vast/scratch/users/jackson.v/retThickness/GWAS
# dataDir=/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/cleanedGeneticFiles/cleanedEURData
  dataDir=/vast/scratch/users/jackson.v/retThickness/GWAS/geneticData/cleanedEURData

cd $workDir

# mkdir -p $workDir/plink
# mkdir -p $workDir/data
# cd $workDir/plink
# wget "https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_x86_64_20221024.zip"
# unzip plink2_linux_x86_64_20221024.zip

mkdir -p /vast/scratch/users/jackson.v/retThickness/GWAS/annot
mkdir -p $workDir/annot/data/
mkdir -p $workDir/annot/fpcData/

cd /vast/scratch/users/jackson.v/retThickness/GWAS/annot
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/chr*sentinels_clumpThresh0.001_withOverlap.txt data
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/sentinels/chr*sentinels_clumpThresh0.001_withOverlap.txt fpcData

rm bonferroniSents.txt
for chr in {1..22} X
do

 awk -F'\t' '$7 < 1.721704e-12 {print $2}' ./data/chr${chr}sentinels_clumpThresh0.001_withOverlap.txt >> bonferroniSents.txt
 awk -F'\t' '{print $2}' ./data/chr${chr}sentinels_clumpThresh0.001_withOverlap.txt >>allSents.txt

done


rm pfilesList*.txt
for chr in {1..22} X
do

../plink/plink2 \
  --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr${chr} \
  --extract allSents.txt \
  --make-pgen \
  --out chr${chr}sentinelsAll 

   echo "chr${chr}sentinelsAll" >> pfilesListAll.txt

   ../plink/plink2 \
  --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr${chr} \
  --extract bonferroniSents.txt \
  --make-pgen \
  --out chr${chr}sentinelsBonf 

   echo "chr${chr}sentinelsBonf" >> pfilesListBonf.txt

 ../plink/plink2 \
   --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr${chr} \
   --extract allSents.txt  \
   --freq \
   --out chr${chr}sentinelsFreq

done

../plink/plink2 \
  --pmerge-list pfilesListAll.txt pfile \
  --make-just-pvar cols=vcfheader \
  --out allChrSentinelsAll

../plink/plink2 \
  --pmerge-list pfilesListBonf.txt pfile \
  --make-just-pvar cols=vcfheader \
  --out allChrSentinelsBonf

rsync -av allChrSentinels*.pvar /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output



## FPC results

rm ./fpc*Sents.txt
for chr in {1..22} X
do

if [ -e "./fpcData/chr${chr}sentinels_clumpThresh0.001_withOverlap.txt" ]; then
 awk -F'\t' '$14 < 5E-8/6 {print $3}' ./fpcData/chr${chr}sentinels_clumpThresh0.001_withOverlap.txt >> fpcBonferroniSents.txt
 awk -F'\t' '{print $3}' ./fpcData/chr${chr}sentinels_clumpThresh0.001_withOverlap.txt >> fpcAllSents.txt

fi

done

rm pfilesListFPC*.txt
for chr in {1..22} X
do

../plink/plink2 \
  --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr${chr} \
  --extract fpcBonferroniSents.txt \
  --make-pgen \
  --out FPCchr${chr}sentinelsBonf 

if [ -e "FPCchr${chr}sentinelsBonf.pvar" ]; then
   echo "FPCchr${chr}sentinelsBonf" >> pfilesListFPCBonf.txt

fi

../plink/plink2 \
  --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr${chr} \
  --extract fpcAllSents.txt \
  --make-pgen \
  --out FPCchr${chr}sentinelsAll 

if [ -e "FPCchr${chr}sentinelsAll.pvar" ]; then
   echo "FPCchr${chr}sentinelsAll" >> pfilesListFPCAll.txt

fi


done

../plink/plink2 \
  --pmerge-list pfilesListFPCBonf.txt pfile \
  --make-just-pvar cols=vcfheader \
  --out FPCallChrSentinelsBonf

../plink/plink2 \
  --pmerge-list pfilesListFPCAll.txt pfile \
  --make-just-pvar cols=vcfheader \
  --out FPCallChrSentinelsAll


rsync -av FPCallChrSentinels*.pvar /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output













# ../plink/plink2 \
#   --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr${chr} \
#   --extract /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/sentinels/chr${chr}sentinelsIDonly_clumpThresh0.001_withOverlap.txt \
#   --make-just-pvar cols=vcfheader \
#   --out chr${chr}FPCsentinels


# ../plink/plink2 \
#   --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr${chr} \
#   --extract /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/chr${chr}sentinelsIDonly_clumpThresh0.001_withOverlap.txt \
#   --freq \
#   --out chr${chr}sentinelsFreq

# ## all SNPs
# ../plink/plink2 \
#   --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr${chr} \
#   --make-just-pvar cols=vcfheader \
#   --out chr${chr}

# ../plink/plink2 \
#   --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr${chr} \
#   --freq \
#   --out chr${chr}Freq

# ## annotate using VEP and concatenate

# # $VEP_PATH/vep \
# #   --input_file  chr${chr}sentinels.pvar \
# #   --format vcf \
# #   --dir_cache ./VEP/ \
# #   --assembly GRCh37 \
# #   --port 3337 \
# #   --cache \
# #   --sift b \
# #   --polyphen b \
# #   --force_overwrite \
# #   --canonical \
# #   --symbol \
# #   --pick \
# #   --gene_phenotype \
# #   --nearest symbol \
# #   --pubmed \
# #   --output_file chr${chr}sentinels_temp.tsv \
# #   --stats_text

# # sed 's/#Uploaded_variation/Uploaded_variation/g' chr${chr}sentinels_temp.tsv | sed '/^##/ d' > chr${chr}sentinels_annotated.tsv
# # rm chr${chr}sentinels_temp.tsv


# ./VEP/ensembl-vep/vep \
#   --input_file  chr${chr}sentinels.pvar \
#   --format vcf \
#   --assembly GRCh37 \
#   --port 3337 \
#   --sift b \
#   --polyphen b \
#   --force_overwrite \
#   --canonical \
#   --symbol \
#   --pick \
#   --gene_phenotype \
#   --nearest symbol \
#   --pubmed \
#   --output_file chr${chr}sentinels_temp.tsv \
#   --stats_text
