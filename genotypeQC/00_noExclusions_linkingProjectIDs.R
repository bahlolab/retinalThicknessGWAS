## Run using R/4.1.2

## Link mactel project and healthy dev aging project IDs -
## cleaned genetic data has latter ids

library(data.table)
library(magrittr)
library(dplyr)


dataDir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/genotypeQCnoExclusions/rawData/"
outDir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/genotypeQCnoExclusions/output/"

mtIDs <- paste0(dataDir,"ukb28541_imp_chr1_v3_s487276.sample") %>%
  fread
setnames(mtIDs, paste0(names(mtIDs),"Mt"))

hdaIDs <- paste0(dataDir,"ukb36610_imp_chr1_v3_s487296.sample") %>%
  fread
setnames(hdaIDs, paste0(names(hdaIDs),"Hda"))

merged <- cbind(mtIDs, hdaIDs)

## sense check
nrow(merged[sexMt==sexHda]) == nrow(merged)

## check mismatches
merged[sexMt!=sexHda]
## withdrawn individuals
##otherwise ok

## read in individuals with cleaned phenotype data
pheno <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fPCscores_noExclusions.csv", select = 1)

mergedPheno <- merged[ID_1Mt %in% pheno[,patID]]

## sense check again
nrow(mergedPheno[sexMt==sexHda]) == nrow(mergedPheno)


## output file with the two IDs

linkFile <- mergedPheno[, .(ID_1Mt, ID_1Hda)] %>%
setnames(., c("patID", "patIDhda"))

fwrite(linkFile, file = paste0(outDir,"idLinkage.txt"), sep ="\t")
