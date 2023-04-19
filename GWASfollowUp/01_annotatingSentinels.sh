#!/bin/bash

#SBATCH -J extractGWsig
#SBATCH -o /vast/scratch/users/jackson.v/retThickness/GWAS/logs/annotatingResults_%A_%a.log
#SBATCH -t 12:0:0
#SBATCH --mem=10GB
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-22


module load htslib
module load perl

export VEP_PATH=/wehisan/bioinf/lab_bahlo/users/jackson.v/resources/VEP/ensembl-vep
export PERL5LIB=$PERL5LIB:${VEP_PATH}/Plugins/
export PATH=~/vcftools/usr/local/bin/:$PATH

chr=$SLURM_ARRAY_TASK_ID

workDir=/vast/scratch/users/jackson.v/retThickness/GWAS
dataDir=/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/cleanedGeneticFiles/cleanedEURData

cd $workDir

# mkdir -p $workDir/plink
# cd $workDir/plink
# wget "https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_x86_64_20221024.zip"
# unzip plink2_linux_x86_64_20221024.zip

# mkdir -p /vast/scratch/users/jackson.v/retThickness/GWAS/annot

cd /vast/scratch/users/jackson.v/retThickness/GWAS/annot

../plink/plink2 \
  --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr${chr} \
  --extract /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/chr${chr}sentinelsIDonly_clumpThresh0.001_withOverlap.txt \
  --make-just-pvar cols=vcfheader \
  --out chr${chr}sentinels


../plink/plink2 \
  --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr${chr} \
  --extract /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/sentinels/chr${chr}sentinelsIDonly_clumpThresh0.001_withOverlap.txt \
  --make-just-pvar cols=vcfheader \
  --out chr${chr}FPCsentinels


../plink/plink2 \
  --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr${chr} \
  --extract /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/chr${chr}sentinelsIDonly_clumpThresh0.001_withOverlap.txt \
  --freq \
  --out chr${chr}sentinelsFreq

## all SNPs
../plink/plink2 \
  --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr${chr} \
  --make-just-pvar cols=vcfheader \
  --out chr${chr}

../plink/plink2 \
  --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr${chr} \
  --freq \
  --out chr${chr}Freq

## annotate using VEP and concatenate

# $VEP_PATH/vep \
#   --input_file  chr${chr}sentinels.pvar \
#   --format vcf \
#   --dir_cache ./VEP/ \
#   --assembly GRCh37 \
#   --port 3337 \
#   --cache \
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

# sed 's/#Uploaded_variation/Uploaded_variation/g' chr${chr}sentinels_temp.tsv | sed '/^##/ d' > chr${chr}sentinels_annotated.tsv
# rm chr${chr}sentinels_temp.tsv


./VEP/ensembl-vep/vep \
  --input_file  chr${chr}sentinels.pvar \
  --format vcf \
  --assembly GRCh37 \
  --port 3337 \
  --sift b \
  --polyphen b \
  --force_overwrite \
  --canonical \
  --symbol \
  --pick \
  --gene_phenotype \
  --nearest symbol \
  --pubmed \
  --output_file chr${chr}sentinels_temp.tsv \
  --stats_text
