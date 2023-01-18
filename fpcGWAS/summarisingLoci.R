#!/usr/bin/env Rscript

library(data.table)
library(magrittr)
library(tidyverse)


loci <- lapply(c(1:22), function(chr) {

results <- lapply(c(1:100), function(fpc) {
  
  print(fpc)
  # slice <- 64
  file <- paste0("/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/fpcGWAS/results/chr",chr,"/chr",chr,"EUR.fpc",fpc,".glm.linear")
  
  fpcResults <- fread(file, select = c("#CHROM", "ID", "POS", "A1_CT",  "BETA", "P"), fill=TRUE) %>%
    setnames(., "#CHROM", "CHR") %>%
    .[P < 5E-8] %>%
    .[, FPC := fpc]
  
  return(fpcResults)
  
}) %>%
  rbindlist


dir <- "/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/fpcGWAS/ClumpedResults"

clumpedFull <- lapply(c(1:100), function(fpc) {
  
  file <- paste0(dir,"/chr",chr,"/chr",chr,"EUR.",fpc,".clumped")

  if(file.exists(file)) {
    fpcResults <- fread(file, fill=T) %>%
      .[, FPC := fpc]
    
    if(nrow(fpcResults) > 0) {
      return(fpcResults)
    }
    
  }
  
}) %>%
  rbindlist %>%
  na.omit




results.filt <- results
clumped <- clumpedFull
sentinels <- NULL
allocatedCheck <- NULL

while(nrow(results.filt)>0) {
  # select top hit
  hit <- results.filt[P==min(results.filt[,P])]
  
  print(paste(hit[,ID], " with P =", hit[,P]))
  
  snp <- hit[1,ID]
  fpc <- hit[1,FPC]
  clumpedResult <- clumped[SNP==snp & FPC==fpc]
  
  if(nrow(clumpedResult)==0) {
    print(paste("SNP", snp, "Not sentinel for chunk"))
    cat("\n")
    
    allFPC <- results.filt[ID==snp, FPC] %>% unique
    check <- data.table(snp=snp,
                        fpc= paste(allFPC, collapse = ","),
                        nFPC = length(allFPC))
    
    allocatedCheck <- rbind(allocatedCheck, check)
    results.filt <- results.filt[!(ID %in% snp & FPC %in% allFPC)]
    
  } else {
    
    tmp <- clumpedResult[,SP2] %>%
      strsplit(., ",") %>%
      unlist
    clumpSNPs <- sapply(strsplit(tmp,"\\("), `[`, 1)
    
    clumpFPCs <- results.filt[ID %in% c(snp,clumpSNPs),FPC] %>% unique
    minFPC <- clumpFPCs %>% min
    n <- length(clumpSNPs)
    noSNPs <- clumpSNPs[1]=="NONE"
    nSNPs <- ifelse(noSNPs, 0, n)
    nFPC <- length(clumpFPCs)
    
    hit <- hit[, nSNPsLocus := nSNPs] %>%
      .[, nFPCLocus := nFPC] %>%
      .[, SNPsInLocus := paste(clumpSNPs, collapse = ",")] %>%
      .[, FPCs := paste(clumpFPCs, collapse = ",")] %>%
      .[, firstFPC := minFPC]
    
    # remove SNPs in clump
    results.filt <- results.filt[!(ID %in% c(snp,clumpSNPs) & FPC %in% c(fpc,clumpFPCs))]
    
    # Add hit onto list of sentinel SNPs
    sentinels <- rbind(sentinels, hit[1,])
    
    print(paste(nSNPs, "additional SNPs in locus"))
    print(paste(nFPC, "additional FPCs significant at locus"))
    print(paste(nrow(results.filt), "SNPs remaining..."))
    cat("\n")
  }
}

return(sentinels)
}) %>%
  rbindlist

fpcsVec <- strsplit(loci[, FPCs], ",") %>% unlist %>% as.numeric()
fpcLoci <- data.table(fpcs = fpcsVec)

ggplot(fpcLoci, aes(x = fpcs)) + geom_bar()
