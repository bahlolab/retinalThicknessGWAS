#!/usr/bin/env Rscript

## Run using R/4.1.2

library(data.table)
library(magrittr)
library(tidyverse)
library(ieugwasr)

sentinels <- lapply(c(1:22, "X"), function(chr) {
  
  chrSent <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/chr",chr,"sentinels_clumpThresh0.001_withOverlap.txt") %>%
    fread() %>%
  .[nSNPsLocus >= 5]
  return(chrSent)
}) %>%
  rbindlist(., idcol = "CHR")

plinkPath <- "/vast/scratch/users/jackson.v/retThickness/GWAS/plink/plink2"

loci <- lapply(c(1:nrow(sentinels)), function(i) {

    locusID <- rep(i, times = sentinels[i, nSNPsLocus]+1)
    sentinelSNPID <- rep(sentinels[i, ID], times = sentinels[i, nSNPsLocus]+1)
    locusSNPID <- sentinels[i, SNPsInLocus] %>% 
        str_split(., ",") %>%  
        unlist %>% 
        c(sentinels[i, ID], .)


    bfilePath <- ifelse(sentinels[i,CHR] == 23,
    "/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/cleanedGeneticFiles/cleanedEURData/plinkBin/EUR_minMaf0.005_minInfo0.8_chrX", 
    paste0("/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/cleanedGeneticFiles/cleanedEURData/plinkBin/EUR_minMaf0.005_minInfo0.8_chr",sentinels[i,CHR]))

     ld <- ld_matrix_local(locusSNPID, 
        bfile = bfilePath, 
        plink_bin = "/stornext/System/data/apps/plink/plink-1.9/bin/plink")

      ldSent <- ld[,colnames(ld) %like% sentinelSNPID[1]]

      ldSNPs <- sapply(strsplit(names(ldSent), "_", fixed = TRUE),
       function(i) paste(head(i, -2), collapse = "_"))

      ldSentDT <- data.table(locusSNPID = ldSNPs, 
                              r2withSentinel = ldSent^2)

    locusDT <- data.table(locusID, sentinelSNPID, locusSNPID) %>%
      .[ldSentDT,  on = "locusSNPID"]   

     return(locusDT) 

}) %>%
rbindlist


fwrite(loci, "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/finalResultsEUR/gwSigLociSNPs.csv", sep =",")

sentinelsOut <- data.table(locusID = 1:nrow(sentinels)) %>%
cbind(., sentinels[, c(1:6,8,9,11)]) %>%
  .[, BonferroniSig := ifelse(P < 5E-8/29041, "Y", "N")]

fwrite(sentinelsOut, "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/finalResultsEUR/gwSigLociSummary.csv", sep =",")
