#!/usr/bin/env Rscript

library(data.table)
library(magrittr)


## Read in list of sentinels meeting genome-wide sig
sentinels <- lapply(c(1:22, "X"), function(chr) {
  
  chrSent <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/chr",chr,"sentinels_clumpThresh0.001_withOverlap.txt") %>%
    fread() %>%
  .[nSNPsLocus >= 5]
  return(chrSent)
}) %>%
  rbindlist(., idcol = "CHR")


## lists of IDs

## All loci meeting 5E-8
gwSigIDs <- sentinels[,ID]

## Restricted to those meeting the bonferroni corrected threshold (ie those reported in the paper).
bonfSigIDs <- sentinels[P < 5E-8/29041, ID]

## read in results for all pixels
## this file contains some SNPs which were filtered out (< 5 SNPs in locus)
results <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/allSentinelsAllPixelsResults_clumpThresh0.001_withOverlap.csv")

## To generate a matrix of betas for all SNPs meeting Bonferroni corrected sig threshold:
mat <- results[ID %in% bonfSigIDs] %>%
  dcast(., pixel ~ ID, value.var = "BETA") %>%
  as.matrix(., rownames = "pixel") 
