#!/bin/bash

cd /vast/scratch/users/jackson.v/retThickness/filteringGeneticData/


mkdir -p plink
cd plink
# wget "https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_x86_64_20220603.zip"
# unzip plink2_linux_x86_64_20220603.zip
wget "https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_x86_64_20221024.zip"
unzip plink2_linux_x86_64_20221024.zip

# project directory
projDir=/vast/scratch/users/jackson.v/retThickness/filteringGeneticData

# specify files

## Withdrawn IDs
myWithdrawnFile=/wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/app36610/sampleLinkage/withdraw36610_137_20221111.txt

## File containing sample QC
mySampleQCFile=/wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/app36610/rawPheno/ukb42082.tab

## SNP QC file
mySNPQCFile=/wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/rawGenetic/QC/ukb_snp_qc.txt

## SNP MAF and info files
mySNPInfoFiles=/wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/rawGenetic/alleleFreqs/ukb_mfi_chr*_v3.txt

## Scripts directory
myScriptsDir=/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/retinalThicknessGWAS/genotypeQC

## directory in lab storage are to copy intermediate files back to
myQCFilesDir=/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/genotypeQC/output

##qctool path
qctoolPath=/wehisan/bioinf/lab_bahlo/users/jackson.v/resources/QCTOOL/qctool_v2.0.1-CentOS6.8-x86_64/


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

module load R/4.1.2
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
mkdir -p ${projDir}/cleanedEURData/
mkdir -p ${projDir}/cleanedCSAData/
mkdir -p ${projDir}/cleanedAFRData/

mkdir -p ${projDir}/cleanedEURData/bgenFilt
mkdir -p ${projDir}/cleanedCSAData/bgenFilt
mkdir -p ${projDir}/cleanedAFRData/bgenFilt

mkdir -p ${projDir}/cleanedEURData/plinkBin
mkdir -p ${projDir}/cleanedCSAData/plinkBin
mkdir -p ${projDir}/cleanedAFRData/plinkBin

mkdir -p ${projDir}/cleanedEURData/plink2Bin
mkdir -p ${projDir}/cleanedCSAData/plink2Bin
mkdir -p ${projDir}/cleanedAFRData/plink2Bin


# filtering of imputed data

## link to all data
for CHR in {1..22};
do
  ln -s /wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/rawGenetic/imputed/ukb_imp_chr${CHR}_v3.bgen data/
done
ln -s /wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/rawGenetic/imputed/ukb_imp_chrX_v3.bgen data/

ln -s /wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/app36610/sampleLinkage/ukb36610_imp_chr1_v3_s487296.sample data/project.sample
ln -s /wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/app36610/sampleLinkage/ukb36610_imp_chrX_v3_s486645.sample data/project_chrX.sample


## run for EURO samples
${myScriptsDir}/filterGeneticData.sh \
 -q  ${qctoolPath} \
 -p  ${projDir}/plink/ \
 -g ./data \
 -s ./data/project.sample \
 -x ./data/project_chrX.sample \
 -k ${myQCFilesDir}/sampleList_doubleIDs_EUR.txt \
 -e ${myQCFilesDir}/snpInclude_minMaf0.005_minInfo0.8.txt \
 -o ${projDir}/cleanedEURData/ \
 -n EUR_minMaf0.005_minInfo0.8

## run for ÇSA samples
${myScriptsDir}/filterGeneticData.sh \
 -q  ${qctoolPath} \
 -p  ${projDir}/plink/ \
 -g ./data \
 -s ./data/project.sample \
 -x ./data/project_chrX.sample \
 -k ${myQCFilesDir}/sampleList_doubleIDs_CSA.txt \
 -e ${myQCFilesDir}/snpInclude_minMaf0.005_minInfo0.8.txt \
 -o ${projDir}/cleanedCSAData/ \
 -n CSA_minMaf0.005_minInfo0.8

## run for AFR samples
${myScriptsDir}/filterGeneticData.sh \
 -q  ${qctoolPath} \
 -p  ${projDir}/plink/ \
 -g ./data \
 -s ./data/project.sample \
 -x ./data/project_chrX.sample \
 -k ${myQCFilesDir}/sampleList_doubleIDs_AFR.txt \
 -e ${myQCFilesDir}/snpInclude_minMaf0.005_minInfo0.8.txt \
 -o ${projDir}/cleanedAFRData/ \
 -n AFR_minMaf0.005_minInfo0.8

## for troubleshooting
qctoolPath= ${qctoolPath}
plinkPath=${projDir}/plink
genDataDir=${projDir}/data
sampFile=${projDir}/data/project.sample
xSampFile=${projDir}/data/project_chrX.sample
keepSamps=${myQCFilesDir}/sampleList_doubleIDs_EUR.txt
extractSNPs=${myQCFilesDir}/snpInclude_minMaf0.005_minInfo0.8.txt
outputDir=${projDir}/cleanedEURData
outName=EUR_minMaf0.005_minInfo0.8



## sense check!

for fi in ./cleaningTemp/plinkErrors/*
do
  echo $fi
  tail -n 1 $fi
done

## copy files to vast data
rsync -av ./cleanedEURData/ /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/cleanedGeneticFiles/cleanedWhiteBritUnrelatedData/
rsync -av ./cleanedCSAdData/ /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/cleanedGeneticFiles/cleanedWhiteBritRelatedData/
rsync -av ./cleanedAFRData/ /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/cleanedGeneticFiles/cleanednonWhiteBritData/

## remove temp files
rm ./cleaningTemp/*
rmdir cleaningTemp
