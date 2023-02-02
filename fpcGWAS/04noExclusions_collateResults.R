#!/usr/bin/env Rscript

library(data.table)
library(magrittr)
library(tidyverse)
library(corrplot)


loci <- lapply(c(1:22), function(chr) {

print(paste("chromosome",chr))

results <- lapply(c(1:10), function(fpc) {
  
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

clumpedFull <- lapply(c(1:10), function(fpc) {
  
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
    
      clumpedResult2 <- clumped[SP2 %like% snp & FPC==fpc]
      clumpSentinel <- clumpedResult2[,SNP]

     if(nrow(sentinels[SNPsInLocus %like% clumpSentinel & FPC %like% fpc])>0) {

      print(paste("SNP", snp, "already assigned to locus"))
      cat("\n")

      check <- data.table(snp=snp,
        POS=clumpedResult2[,BP],
        clumpSentinel=clumpSentinel,
        FPCs=paste(all.fpcs, collapse = ","),
        allocated = "Y",
        otherClumpSNPs=clumpedResult2[,SP2])

      allocatedCheck <- rbind(allocatedCheck, check)


     } else {

      check <- data.table(snp=snp,
        POS=clumpedResult2[,BP],
        clumpSentinel=clumpSentinel,
        FPCs=paste(all.fpcs, collapse = ","),
        allocated = "N",
        otherClumpSNPs=clumpedResult2[,SP2])
      allocatedCheck <- rbind(allocatedCheck, check)

     }

    results.filt <- results.filt[!(ID %in% snp & FPC %like% fpc)]

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
    print(paste(nrow(results.filt), "SNPs remaining..."))
    cat("\n")
}
}


fwrite(sentinels, file =paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/sentinels/chr",chr,"sentinels.txt"), sep = "\t")
fwrite(sentinels[,.(ID)], , file =paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/sentinels/chr",chr,"sentinelsIDonly.txt"), col.names=F)





mat <- results[ID %in% sentinels[,ID]] %>%
dcast(., FPC ~ ID, value.var = "BETA", fun.aggregate=mean) %>%
  as.matrix(., rownames = "FPC")
mat[is.na(mat)] <- 0
class(mat) <- "numeric"

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/sentinels/chr",chr,"_senitnelsBetasCorrPlot.png"), width = 1200, height = 1200)
mat %>% cor %>% corrplot
dev.off()

return(sentinels)
}) %>%
rbindlist(., idcol = "chr")

fwrite(loci, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/sentinels/allChr_sentinels.txt")