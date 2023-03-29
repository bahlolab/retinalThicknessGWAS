library(data.table)
library(magrittr)
library(tidyverse)
library(corrplot)
library(patchwork)
library(biomaRt)


fpcSentinels <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/sentinels/allChr_sentinel_clumpThresh0.001_withOverlap.txt")

sentinels <- lapply(c(1:22), function(chr) {
  
  chrSent <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/chr",chr,"sentinels_clumpThresh0.001_withOverlap.txt") %>%
    fread()
  return(chrSent)
}) %>%
  rbindlist(., idcol = "CHR")

fpcSentinelsAllSNPs <- c(fpcSentinels[,ID], fpcSentinels[,SNPsInLocus]) %>% 
  strsplit(., ",") %>%
  unlist

sentinelsAllSNPs <- c(sentinels[,ID], sentinels[,SNPsInLocus]) %>% 
  strsplit(., ",") %>%
  unlist

fpcSentinels %>% nrow()
fpcSentinels[ID %in% sentinelsAllSNPs] %>% nrow()

sentinels %>% nrow()
sentinels[ID %in% fpcSentinelsAllSNPs] %>% nrow()


## plot FPC results which aren't amongst pixelwise hits.

fpcMismatch <- fpcSentinels[!ID %in% sentinelsAllSNPs] %>%
  .[, chr := NULL]

check <- lapply(c(1:22), function(chr) {
  
  dir <- paste0("/vast/scratch/users/jackson.v/retThickness/GWAS/GWsigResults/chr",chr)
  
  print(paste("chromosome",chr))
  
  slice <- 1
  file <- paste0(dir,"/chr",chr,"Slice",slice,"_5e-5Sig.txt")
  
    slice1 <- fread(file) %>%
      setnames(., "#POS", "POS") %>%
      .[ID %in% fpcMismatch[,ID]]
  
  if(nrow(slice1 > 0)) {
    
  chrResults <- lapply(c(2:119), function(slice) {
    
    #print(paste(slice))
    
    file <- paste0(dir,"/chr",chr,"Slice",slice,"_5e-5Sig.txt")
    
    if(file.exists(file)) {
      sliceResults <- fread(file) %>%
        setnames(., "#POS", "POS") %>%
        .[ID %in% fpcMismatch[,ID]]
      return(sliceResults)
      
    } else {
      print(paste("slice",slice,"missing"))
    }
    
  }) %>%
    rbindlist %>%
    rbind(slice1, .) %>%
    .[, chrom := chr] %>%
    .[, log10P := (-1)*log10(P)] %>%
    .[, c("y", "x") := tstrsplit(pixel, "_", type.convert=TRUE)]
  
  
  sentinels <- fpcMismatch[CHR %in% chr,ID] %>% unique
  
  lapply(sentinels, function(snp) {
    
    topFPC <- fpcMismatch[ID==snp, FPC]
    otherFPC <- fpcMismatch[ID==snp, clumpFPCs]
    
    
    snpResult <- chrResults[ID == snp]
    
    statsPlots <- lapply(c("BETA", "T_STAT", "log10P"), function(stat) {
      plot <- ggplot(snpResult) +
        geom_tile(aes_string(x = "x", y = "y", fill = stat)) +
        # geom_path(aes(x=col,y=row,color=area),data = areas, size = 0.5) +
        scale_fill_gradient2() +
        scale_y_reverse() +
        theme_bw() +
        theme(legend.position = "bottom")+
        ggtitle(paste(snp,"- topFPC:",topFPC,"; other FPCs:",otherFPC))
      
      return(plot)
    })
    png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/fpcSigPixelNot/",snp,"AssocsPixelwise.png"), width = 1800, height = 600)
    reduce(statsPlots , `+`) %>%
      print
    dev.off()
    
  })
  
  return(chrResults)
  }
}) %>%
  rbindlist   



## annotate pixel wise results with FPC hits, and previous hits.

rsids <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/misc/rsids_Gao.txt") %>%
  .[, start := NULL]
gao <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/misc/macularThicknessLoci_Gao.txt") %>%
rsids[., on = c("chr", "pos")]

currant1 <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/misc/retinalThicknessLoci_Currant.csv")
currant2 <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/misc/retinalThicknessLoci_Currant_2023.txt")

gaoPos <- list(gao[1:2, chr], gao[1:2, pos])

gaoAnno <- getBM(attributes = c('refsnp_id', 'allele', 'chr_name', 'chrom_start'), 
      filters = c('chr_name', 'start'), 
      values = gaoPos, 
      mart = snp_mart)


snps <- sentinels[, ID]

pixelAnnotated <- lapply(snps, function(id) {

result <- sentinels[ID==id]

allLocusSNPs <-  c(result[,ID], result[,SNPsInLocus]) %>% 
  strsplit(., ",") %>%
  unlist

fpcMatch <- fpcSentinels[ID %in% allLocusSNPs]

if(nrow(fpcMatch) > 0) {

  fpcOut <- data.table(fpcSig = "Y"
    topFPC = fpcMatch[,FPC]
    allFPCs = fpcMatch[, clumpFPCs]
    sameSentinel = ifelse(fpcMatch[,ID]==id, "Y", "N")
    FPCsentinel = fpcMatch[,ID]
    FPCtopP = fpcMatch[,P])

} else {
  fpcOut <- data.table(fpcSig = "N"
  topFPC = NA
  allFPCs = NA
  sameSentinel = NA
  FPCsentinel = NA
  FPCtopP = NA)

}



})


  
topSents <- highConf[,ID]

highConfFPC <- lapply(topSents, function(id) {
  
  fpcResult <- fpcSentinels[ID==id | SNPsInLocus %like% id] %>%
    .[, sameSentinel := ifelse(ID==id, "Y", "N")] %>%
    .[, .(ID, REF, ALT, A1_FREQ, BETA, P, FPC,clumpFPCs, sameSentinel)] %>%
    setnames(., "ID", "fpcSentinel") %>%
    setnames(., "BETA", "fpcBETA") %>%
    setnames(., "P", "fpcP") %>%
    setnames(., "FPC", "topFPC") 
  
  outResult <- cbind(highConf[ID==id], fpcResult) %>%
    .[, .(CHR, POS, ID, REF, ALT, A1, A1_FREQ, BETA, SE, T_STAT, P, pixel, nSNPsLocus, nPixelsLocus, fpcSentinel, fpcBETA, fpcP, topFPC, sameSentinel, clumpFPCs, SNPsInLocus, pixels)]
    
    
})  %>%
  rbindlist

highConfFPC[order(P), 1:20] %>% head(., n=10)

highConfFPC[order(fpcP), 1:20] %>% head(., n=10)





























check <- check %>%
  .[, log10P := (-1)*log10(P)] %>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert=TRUE)]


sentinels <- fpcMismatch[,ID] %>% unique

lapply(sentinels, function(snp) {
  
  topFPC <- fpcMismatch[ID==snp, fpc]
  otherFPC <- fpcMismatch[ID==snp, fp]
  
  
  snpResult <- check[ID == snp]
  
  statsPlots <- lapply(c("BETA", "T_STAT", "log10P"), function(stat) {
    plot <- ggplot(snpResult) +
      geom_tile(aes_string(x = "x", y = "y", fill = stat)) +
      # geom_path(aes(x=col,y=row,color=area),data = areas, size = 0.5) +
      scale_fill_gradient2() +
      scale_y_reverse() +
      theme_bw() +
      theme(legend.position = "bottom")+
      ggtitle(paste(snp,"- topFPC:"topFPC,"; other FPCs:",otherFPC))
    
    return(plot)
  })
  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/fpcMismatches/",snp,"AssocsPixelwise.png"), width = 1800, height = 600)
  reduce(statsPlots , `+`) %>%
    print
  dev.off()
  
})
