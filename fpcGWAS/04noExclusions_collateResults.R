#!/usr/bin/env Rscript

library(data.table)
library(magrittr)
library(tidyverse)
library(corrplot)


loci <- lapply(c(1:22, "X"), function(chr) {

print(paste("chromosome",chr))

results <- lapply(c(1:6), function(fpc) {
  
  # print(fpc)
  # slice <- 64
  file <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/results/chr",chr,"/chr",chr,"EUR.fpc",fpc,".glm.linear")
  
  fpcResults <- fread(file) %>%
    setnames(., "#CHROM", "CHR") %>%
    .[, FPC := fpc]
  
  return(fpcResults)
  
}) %>%
  rbindlist


dir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/clumpedResults"

clumpedFull <- lapply(c(1:6), function(fpc) {
  
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


if(nrow(clumpedFull ) > 0) {

results.filt <- results[P<5e-5]
clumped <- clumpedFull
sentinels <- NULL
allocatedCheck <- NULL



while(nrow(results.filt[P<5e-8])>0) {
		# select top hit
		hit <- results.filt[P==min(results.filt[,P])]

    print(paste(hit[,ID], " with P =", hit[,P]))

    snp <- hit[1,ID]
    fpc <- hit[1,FPC]
    clumpedResult <- clumped[SNP==snp & FPC==fpc]

    all.fpcs <- results.filt[ID==snp, FPC] %>% unique

    if(nrow(clumpedResult)==0) {
      print(paste("SNP", snp, "Not sentinel for chunk"))
      cat("\n")
    
       # find clump
      fpcClumps  <- clumped[FPC==fpc]

      whichClump <- lapply(fpcClumps[,SNP], function(sent) {

          tmp <- fpcClumps[SNP==sent,SP2] %>%
            strsplit(., ",") %>%
            unlist
          clumpSNPs <- sapply(strsplit(tmp,"\\("), `[`, 1)

        inClump <-  ifelse(snp %in% clumpSNPs, T, F)
        return(inClump)

      }) %>%
      unlist 

      clumpedResult2 <- fpcClumps[whichClump]
      clumpSentinel <- clumpedResult2[,SNP]
       tmp <- fpcClumps[whichClump,SP2] %>%
          strsplit(., ",") %>%
          unlist
        clumpSNPs <- sapply(strsplit(tmp,"\\("), `[`, 1)


      whichSentinel <- lapply(sentinels[,ID], function(sent) {

          tmp <- sentinels[ID==sent, SNPsInLocus ] %>%
            strsplit(., ",") %>%
            unlist
          sentClumpSNPs <- sapply(strsplit(tmp,"\\("), `[`, 1)

        inClump <-  ifelse(any(sent %in% clumpSentinel | sentClumpSNPs %in% clumpSentinel), T, F)
        return(inClump)

      }) %>%
      unlist 

     if(any(whichSentinel)) {

      print(paste("SNP", clumpSentinel, "is top SNP for this clump"))
      print(paste("Sentinel", sentinels[whichSentinel, ID], "also within this clump"))
      print(paste("removing other clump SNPs;", nrow(results.filt[!(ID %in% clumpSNPs & FPC %in% all.fpcs) & P<5E-8]), "SNPs remaining..."))
      cat("\n")

      check <- data.table(snp=snp,
        clumpTopSNPPOS=clumpedResult2[,BP],
        clumptopSNP=clumpSentinel,
        topFPC=fpc,
        nFPC = length(all.fpcs),
        assignedSentinelPOS = sentinels[whichSentinel, POS],
        assignedSentinelSNP = sentinels[whichSentinel, ID],
        assignedSentinelFPC = sentinels[whichSentinel, FPC],
        allocated = "Y",
        FPCs=paste(all.fpcs, collapse = ","),
        otherClumpSNPs=clumpedResult2[,SP2])

      allocatedCheck <- rbind(allocatedCheck, check)


     } else {

      check <- data.table(snp=snp,
        clumpTopSNPPOS=clumpedResult2[,BP],
        clumptopSNP=clumpSentinel,
        topFPC=fpc,
        nFPC = length(all.fpcs),
        assignedSentinelPOS = NA,
        assignedSentinelSNP = NA,
        assignedSentinelFPC = NA,
        allocated = "N",
        FPCs=paste(all.fpcs, collapse = ","),
        otherClumpSNPs=clumpedResult2[,SP2])

      allocatedCheck <- rbind(allocatedCheck, check)

     }

    results.filt <- results.filt[!(ID %in% clumpSNPs & FPC %in% all.fpcs)]

    } else {

    tmp <- clumpedResult[,SP2] %>%
    strsplit(., ",") %>%
    unlist
    clumpSNPs <- sapply(strsplit(tmp,"\\("), `[`, 1)

    clumpFPCs <- results.filt[ID %in% c(snp,clumpSNPs), FPC] %>% unique
    n <- length(clumpSNPs)
    noSNPs <- clumpSNPs[1]=="NONE"
    nSNPs <- ifelse(noSNPs, 0, n)
    nFPCs <- length(clumpFPCs)

		hit <- hit[, nSNPsLocus := nSNPs] %>%
      .[, nFPCsLocus := nFPCs] %>%
      .[, SNPsInLocus := paste(clumpSNPs, collapse = ",")] %>%
     .[, clumpFPCs := paste(clumpFPCs, collapse = ",")]

		# remove SNPs in clump
		results.filt <- results.filt[!(ID %in% c(snp,clumpSNPs) & FPC %in% c(fpc,clumpFPCs))]

		# Add hit onto list of sentinel SNPs
		sentinels <- rbind(sentinels, hit[1,])

    print(paste(nSNPs, "additional SNPs in locus"))
    print(paste(nFPCs, "additional FPCs significant at locus"))
    print(paste(nrow(results.filt[P<5E-8]), "SNPs remaining..."))
    cat("\n")
}
}


fwrite(sentinels, file =paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/sentinels/chr",chr,"sentinels_clumpThresh0.001_withOverlap.txt"), sep = "\t")
fwrite(sentinels[,.(ID)], file =paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/sentinels/chr",chr,"sentinelsIDonly_clumpThresh0.001_withOverlap.txt"), col.names=F)

if(!is.null(allocatedCheck)) {
  fwrite(allocatedCheck, file = paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/sentinels/chr",chr,"_checkAllocated_clumpThresh0.001_withOverlap.txt"), sep = "\t")
}




mat <- results[ID %in% sentinels[,ID]] %>%
dcast(., FPC ~ ID, value.var = "BETA", fun.aggregate=mean) %>%
  as.matrix(., rownames = "FPC")
mat[is.na(mat)] <- 0
class(mat) <- "numeric"

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/sentinels/chr",chr,"_senitnelsBetasCorrPlot.png"), width = 1200, height = 1200)
mat %>% cor %>% corrplot
dev.off()

return(sentinels)
}
}) %>%
rbindlist(., idcol = "chr")

fwrite(loci, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/sentinels/allChr_sentinel_clumpThresh0.001_withOverlap.txt")

