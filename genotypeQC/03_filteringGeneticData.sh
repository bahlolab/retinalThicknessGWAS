#!/bin/bash

cd /vast/scratch/users/jackson.v/retThickness/filteringGeneticData/


## download software
# get it
wget http://code.enkre.net/bgen/tarball/release/bgen.tgz
tar -xvf bgen.tgz
mv bgen.tgz/ bgen
cd bgen
# compile it
./waf configure
./waf
# test it
./build/test/unit/test_bgen


## qctool - use version in user area
# wget https://code.enkre.net/qctool/zip/release/qctool.tgz
# unzip qctool.tgz
#
# cd qctool
# ./waf-1.5.18 configure
# ./waf-1.5.18
#
# ./build/release/qctool_v2.0-<version> -help
#
#
# wget https://www.well.ox.ac.uk/~gav/resources/qctool_v2.2.0-osx.tgz
# tar -xvf qctool_v2.2.0-osx.tgz
#

# project directory
projDir=/vast/scratch/users/jackson.v/retThickness/filteringGeneticData

# specify files

## Withdrawn IDs
myWithdrawnFile=/wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/app36610/sampleLinkage/withdraw36610_120_20220525.txt

## File containing sample QC
mySampleQCFile=/wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/app36610/rawPheno/ukb42082.tab

## SNP QC file
mySNPQCFile=/wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/rawGenetic/QC/ukb_snp_qc.txt

## SNP MAF and info files
mySNPInfoFiles=/wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/rawGenetic/alleleFreqs/ukb_mfi_chr*_v3.txt

## Scripts directory
myScriptsDir=/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/genotypeQC/scripts

## directory in lab storage are to copy intermediate files back to
myQCFilesDir=/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/genotypeQC/output

## specify plink path
plinkPath=/wehisan/bioinf/lab_bahlo/users/jackson.v/resources/plink_linux_x86_64_20190617/

########################################################################################################################
## Set up directory structure
mkdir ${projDir}/data
mkdir ${projDir}/qcFiles

# link to data
ln -s ${myWithdrawnFile}  ${projDir}/data/withdrawn.csv
ln -s ${mySampleQCFile}  ${projDir}/data/sampleQC.tab
ln -s ${mySNPQCFile}  ${projDir}/data/
ln -s ${mySNPInfoFiles}  ${projDir}/data/

########################################################################################################################
# Run cleaning R Scripts

module load R/4.0.2
## set directory
cd ${projDir}

## Run sample cleaning scripts
${myScriptsDir}/sampleQC.R \
  --withdrawn data/withdrawn.csv \
  --sampQC data/sampleQC.tab \
  --outDir $myQCFilesDir

## Run SNP cleaning scripts
${myScriptsDir}/variantQC.R \
  --snpQC data/ukb_snp_qc.txt \
  --mfiDir data \
  --outDir $myQCFilesDir
# rsync -av qcFiles $myQCFilesDir

######################################################################################################################
# Filtering and chunking Scripts

## Set up directory structure
mkdir ${projDir}/cleanedWhiteBritUnrelatedData/
mkdir ${projDir}/cleanedWhiteBritUnrelatedData/plink
mkdir ${projDir}/cleanedWhiteBritUnrelatedData/imputed

mkdir ${projDir}/cleanedWhiteBritRelatedData/
mkdir ${projDir}/cleanedWhiteBritRelatedData/plink
mkdir ${projDir}/cleanedWhiteBritRelatedData/imputed

mkdir ${projDir}/cleanednonWhiteBritData/
mkdir ${projDir}/cleanednonWhiteBritData/plink
mkdir ${projDir}/cleanednonWhiteBritData/imputed

# Running Plink filtering to create file for GRM


## link to all data

for CHR in {1..22};
do
  ln -s /wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/rawGenetic/plink/ukb_cal_chr${CHR}_v2.bed data/
  ln -s /wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/rawGenetic/plink/ukb_snp_chr${CHR}_v2.bim data/
done

ln -s /wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/app36610/sampleLinkage/ukb36610_cal_chr1_v2_s488264.fam data/project.fam

## generate list of plink fils to merge
for CHR in {1..22};
do
  echo ${projDir}/data/ukb_cal_chr${CHR}_v2.bed ${projDir}/data/ukb_snp_chr${CHR}_v2.bim ${projDir}/data/project.fam >> data/plinkList.txt
done

## run for related WhiteBrit samples only
${myScriptsDir}/filterPlinkData.sh \
 -p ${plinkPath} \
 -l ./data/plinkList.txt \
 -k ${myQCFilesDir}/sampleIncludeUnrelatedWhiteBrit_plink.txt \
 -e ${myQCFilesDir}/snpsInclude_grm_plink.txt \
 -o ./cleanedWhiteBritUnrelatedData/plink/grmSNPsUnrelatedWhiteBrit

 ${myScriptsDir}/filterPlinkData.sh \
  -p ${plinkPath} \
  -l ./data/plinkList.txt \
  -k ${myQCFilesDir}/sampleIncludeRelatedWhiteBrit_plink.txt \
  -e ${myQCFilesDir}/snpsInclude_grm_plink.txt \
  -o ./cleanedWhiteBritRelatedData/plink/grmSNPsRelatedWhiteBrit

  ${myScriptsDir}/filterPlinkData.sh \
   -p ${plinkPath} \
   -l ./data/plinkList.txt \
   -k ${myQCFilesDir}/sampleIncludeUnrelatedNonWBIDs_plink.txt \
   -e ${myQCFilesDir}/snpsInclude_grm_plink.txt \
   -o ./cleanednonWhiteBritData/plink/grmSNPsUnrelatedNonWB



# filtering of imputed data

## specify qctool path
qctoolPath=/wehisan/bioinf/lab_bahlo/users/jackson.v/resources/QCTOOL/qctool_v2.0.1-CentOS6.8-x86_64/

## link to all data
for CHR in {1..22};
do
  ln -s /wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/rawGenetic/imputed/ukb_imp_chr${CHR}_v3.bgen data/
done
ln -s /wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/rawGenetic/imputed/ukb_imp_chrX_v3.bgen data/

ln -s /wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/app36610/sampleLinkage/ukb36610_imp_chr1_v3_s487296.sample data/project.sample
ln -s /wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/app36610/sampleLinkage/ukb36610_imp_chrX_v3_s486645.sample data/project_chrX.sample


## run for WhiteBrit samples
${myScriptsDir}/filterImputedData.sh \
 -q ${qctoolPath} \
 -g ./data \
 -s ./data/project.sample \
 -x ./data/project_chrX.sample \
 -k ${myQCFilesDir}/sampleIncludeUnrelatedWhiteBrit.txt \
 -e ${myQCFilesDir}/snpIncludeAltID_minMaf0.001_minInfo0.5.txt \
 -o ./cleanedWhiteBritUnrelatedData/imputed/ukbb_minMaf0.001_minInfo0.5 \
 -n UnrelatedWhiteBrit

 ${myScriptsDir}/filterImputedData.sh \
  -q ${qctoolPath} \
  -g ./data \
  -s ./data/project.sample \
  -x ./data/project_chrX.sample \
  -k ${myQCFilesDir}/sampleIncludeRelatedWhiteBrit.txt \
  -e ${myQCFilesDir}/snpIncludeAltID_minMaf0.001_minInfo0.5.txt \
  -o ./cleanedWhiteBritRelatedData/imputed/ukbb_minMaf0.001_minInfo0.5 \
  -n RelatedWhiteBrit

  ${myScriptsDir}/filterImputedData.sh \
   -q ${qctoolPath} \
   -g ./data \
   -s ./data/project.sample \
   -x ./data/project_chrX.sample \
   -k ${myQCFilesDir}/sampleIncludeUnrelatedNonWBIDs.txt \
   -e ${myQCFilesDir}/snpIncludeAltID_minMaf0.001_minInfo0.5.txt \
   -o ./cleanednonWhiteBritData/imputed/ukbb_minMaf0.001_minInfo0.5 \
   -n UnrelatedNonWB

## for troubleshooting
# qctoolPath=${qctoolPath}
# genDataDir=./data
# sampFile=./data/project.sample
# xSampFile=./data/project_chrX.sample
# keepSamps=${myQCFilesDir}/sampleIncludeUnrelatedWhiteBrit.txt
# extractSNPs=${myQCFilesDir}/snpIncludeAltID_minMaf0.001_minInfo0.8.txt
# outputPrefix=./cleanedWhiteBritUnrelatedData/imputed/ukbb_minMaf0.0001_minInfo0.8
# jobName=UnrelatedWhiteBrit

## sense check!

for fi in ./cleaningTemp/qctoolErrors/*
do
  echo $fi
  tail -n 1 $fi
done

## copy files to vast data
rsync -av ./cleanedWhiteBritUnrelatedData/ /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/cleanedGeneticFiles/cleanedWhiteBritUnrelatedData/
rsync -av ./cleanedWhiteBritRelatedData/ /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/cleanedGeneticFiles/cleanedWhiteBritRelatedData/
rsync -av ./cleanednonWhiteBritData/ /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/cleanedGeneticFiles/cleanednonWhiteBritData/

## remove temp files
rm ./cleaningTemp/*
rmdir cleaningTemp
