#!/bin/bash

mkdir -p /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions 
mkdir -p /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/rawData
mkdir -p /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output
mkdir -p /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS
mkdir -p /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/results
mkdir -p /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/ClumpedResults

for chr in {1..22}
do

mkdir -p /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/results/chr${chr}/
mkdir -p /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/ClumpedResults/chr${chr}/

done


  cd /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions

  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/genotypeQCnoExclusions/output/idLinkage.txt rawData
  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fPCscores_noExclusions.csv rawData
  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/processedData/scansUnadjustedFinal.csv rawData
  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/processedData/ rawData
  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/app28541/phenoData/ukb41258.tab rawData
  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/app36610/rawPheno/ukb32825.tab rawData
  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/genotypeQCnoExclusions/output/sampleList_singleIDs_EUR.txt rawData
  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/genotypeQCnoExclusions/output/sampleList_singleIDs_AFR.txt rawData
  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/genotypeQCnoExclusions/output/sampleList_singleIDs_CSA.txt rawData

  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/app36610/rawPheno/Return2442/all_pops_non_eur_pruned_within_pop_pc_covs.tsv rawData
  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/app36610/rawPheno/Return2442/ukb36610bridge31063.txt rawData
  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/genotypeQCnoExclusions/output/sampleGenoQC.csv rawData
