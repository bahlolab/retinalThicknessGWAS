#!/bin/bash

  cd /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWAS

  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/genotypeQC/output/idLinkage.txt rawData
  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/scansUnadjustedFinalfPCexclusions.csv rawData
  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/app28541/phenoData/ukb41258.tab rawData
  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/app36610/rawPheno/ukb32825.tab rawData
  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/app36610/rawPheno/Return2442/all_pops_non_eur_pruned_within_pop_pc_covs.tsv rawData
  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/app36610/rawPheno/Return2442/ukb36610bridge31063.txt rawData
  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/genotypeQC/output/sampleGenoQC.csv rawData
  rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/processedData/fpcPhenotypes.txt rawData
