#!/bin/bash

mkdir -p /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/sensitivityAnalysesGWAS
mkdir -p /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/sensitivityAnalysesGWAS/rawData
mkdir -p /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/sensitivityAnalysesGWAS/output
mkdir -p /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/sensitivityAnalysesGWAS/output/plots


cd /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/sensitivityAnalysesGWAS

rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/genotypeQCnoExclusions/output/idLinkage.txt rawData
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/processedData/scansUnadjustedFinal.csv rawData
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/app28541/phenoData/ukb41258.tab rawData
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/app36610/rawPheno/ukb32825.tab rawData
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/app28541/phenoData/ukb11226.tab rawData
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/genotypeQC/output/sampleGenoQC.csv rawData
rsync -av /stornext/Bioinf/data/lab_bahlo/projects/misc/UKBiobank/data/app28541/phenoData/smokingStatusIDs/* rawData
 
  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/genotypeQCnoExclusions/output/sampleList_singleIDs_EUR.txt rawData
  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/genotypeQCnoExclusions/output/sampleList_singleIDs_AFR.txt rawData
  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/genotypeQCnoExclusions/output/sampleList_singleIDs_CSA.txt rawData

  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/app36610/rawPheno/Return2442/all_pops_non_eur_pruned_within_pop_pc_covs.tsv rawData
  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/app36610/rawPheno/Return2442/ukb36610bridge31063.txt rawData
  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/genotypeQCnoExclusions/output/sampleGenoQC.csv rawData
