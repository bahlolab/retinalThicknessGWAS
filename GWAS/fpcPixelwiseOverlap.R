library(data.table)
library(magrittr)
library(tidyverse)
library(corrplot)
library(patchwork)


fpcSentinels <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/sentinels/allChr_sentinels.txt")

sentinels <- lapply(c(1:22), function(chr) {
  
  chrSent <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/chr",chr,"sentinels_clumpThresh0.001.txt") %>%
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
    png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/fpcMismatches/",snp,"AssocsPixelwise.png"), width = 1800, height = 600)
    reduce(statsPlots , `+`) %>%
      print
    dev.off()
    
  })
  
  return(chrResults)
  }
}) %>%
  rbindlist   




highConf <- sentinels[ID %in% fpcSentinelsAllSNPs] 
  
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

