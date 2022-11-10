#!/bin/bash

  cd /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWAS

  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/genotypeQC/output/idLinkage.txt rawData
  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/processedData/scansUnadjustedFinal.csv rawData
  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/app28541/phenoData/ukb41258.tab rawData
  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/app36610/rawPheno/ukb32825.tab rawData
  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/genotypeQC/output/sampleGenoQC.csv rawData
  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/processedData/fpcPhenotypes.txt rawData
