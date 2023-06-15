#!/bin/bash

cd /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenoAssociations

rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/genotypeQC/output/idLinkage.txt rawData
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/processedData/scansUnadjustedFinal.csv rawData
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/app28541/phenoData/ukb41258.tab rawData
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/genotypeQC/output/sampleGenoQC.csv rawData
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/genotypeQCnoExclusions/output/sampleList_singleIDs_*.txt rawData
