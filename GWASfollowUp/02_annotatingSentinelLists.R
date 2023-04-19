library(data.table)
library(magrittr)
library(tidyverse)
library(corrplot)
library(patchwork)
library(UpSetR)


fpcSentinels <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/sentinels/allChr_sentinel_clumpThresh0.001_withOverlap.txt")

sentinels <- lapply(c(1:22, "X"), function(chr) {
  
  chrSent <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/chr",chr,"sentinels_clumpThresh0.001_withOverlap.txt") %>%
    fread()
  return(chrSent)
}) %>%
  rbindlist(., idcol = "CHR")

fpcSentinelsAllSNPs <- c(fpcSentinels[,ID], fpcSentinels[,SNPsInLocus]) %>% 
  strsplit(., ",") %>%
  unlist

fpcSentinelsAllSNPsBon <- c(fpcSentinels[P < (5E-8/6),ID], fpcSentinels[P < (5E-8/6),SNPsInLocus]) %>% 
  strsplit(., ",") %>%
  unlist

sentinelsAllSNPs <- c(sentinels[,ID], sentinels[P < (5E-8/29041),SNPsInLocus]) %>% 
  strsplit(., ",") %>%
  unlist

sentinelsAllSNPsBon <- c(sentinels[P < (5E-8/29041),ID], sentinels[P < (5E-8/29041),SNPsInLocus]) %>% 
  strsplit(., ",") %>%
  unlist

## all overlap
fpcSentinels %>% nrow()
fpcSentinels[ID %in% sentinelsAllSNPs] %>% nrow()

# Bon corrected overlap
fpcSentinels[P < (5E-8/6)] %>% nrow()
fpcSentinels[P < (5E-8/6) & ID %in% sentinelsAllSNPsBon] %>% nrow()
fpcSentinels[P < (5E-8/6) & !(ID %in% sentinelsAllSNPsBon)] %>% nrow()

## all overlap
sentinels %>% nrow()
sentinels[ID %in% fpcSentinelsAllSNPs] %>% nrow()

# Bon corrected overlap
sentinels[P < (5E-8/29041)] %>% nrow()
sentinels[P < (5E-8/29041) & ID %in% fpcSentinelsAllSNPsBon] %>% nrow()
sentinels[P < (5E-8/29041) & !(ID %in% fpcSentinelsAllSNPsBon)] %>% nrow()


## plot FPC results which aren't amongst pixelwise hits.

fpcMismatch <- fpcSentinels[!ID %in% sentinelsAllSNPs] %>%
  .[, chr := NULL]

check <- lapply(c(1:22, "X"), function(chr) {
  
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


## pixelwise
vep <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/misc/pixelWiseSentinelsVEP20230327.txt") %>%
  .[, c("#Uploaded_variation", "Consequence", "SYMBOL", "Existing_variation", "DISTANCE", "AF", "PUBMED", "PHENOTYPES")] %>%
  setnames(., "#Uploaded_variation", "ID")


freqs <- lapply(c(1:22), function(chr) {
  fread(paste0("/vast/scratch/users/jackson.v/retThickness/GWAS/annot/chr",chr,"sentinelsFreq.afreq"))
}) %>%
  rbindlist %>%
  .[, .(ID, REF, ALT, ALT_FREQS)]

snps <- sentinels[, ID]

pixelsAnnotated <- lapply(snps, function(id) {
  
  result <- sentinels[ID==id]
  
  allLocusSNPs <-  c(result[,ID], result[,SNPsInLocus]) %>% 
    strsplit(., ",") %>%
    unlist
  
  fpcMatch <- fpcSentinels[ID %in% allLocusSNPs]
  
  if(nrow(fpcMatch) > 0) {
    
    fpcOut <- data.table(fpcSig = "Y",
                         topFPC = paste(fpcMatch[,FPC], collapse = ";"),
                         allFPCs = paste(fpcMatch[, clumpFPCs], collapse = ";"),
                         sameSentinel = ifelse(id %in% fpcMatch[,ID], "Y", "N"),
                         FPCsentinel = paste(fpcMatch[,ID], collapse = ";"),
                         FPCtopP = paste(fpcMatch[,P], collapse = ";"))
    
  } else {
    fpcOut <- data.table(fpcSig = "N",
                         topFPC = NA,
                         allFPCs = NA,
                         sameSentinel = NA,
                         FPCsentinel = NA,
                         FPCtopP = NA)
    
  }
  
  gaoMatch <- gao[rsid %in% allLocusSNPs]
  
  if(nrow(gaoMatch) > 0) {
    
    inGao <- data.table(gaoSig = "Y",
                        sameSentinelGao = ifelse(id %in% gaoMatch[,rsid], "Y", "N"))
    
  } else {
    inGao <- data.table(gaoSig = "N",
                        sameSentinelGao = NA)
  }
  
  currantInnerMatch <- currant1[SNP %in% allLocusSNPs]
  
  if(nrow(currantInnerMatch) > 0) {
    
    inCurrantInner <- data.table(currantInnerSig = "Y",
                                 sameSentinelCurrantInner = ifelse(id %in% currantInnerMatch[,SNP], "Y", "N"))
    
  } else {
    inCurrantInner <- data.table(currantInnerSig = "N",
                                 sameSentinelCurrantInner = NA)
  }
  
  currantOuterMatch <- currant2[SNP %in% allLocusSNPs]
  
  if(nrow(currantOuterMatch) > 0) {
    
    inCurrantOuter <- data.table(currantOuterSig = "Y",
                                 sameSentinelCurrantOuter = ifelse(any(currantOuterMatch[,SNP] %in% id), "Y", "N"))
    
  } else {
    inCurrantOuter <- data.table(currantOuterSig = "N",
                                 sameSentinelCurrantOuter = NA)
  }
  
  snpResult <- sentinels[ID==id, 1:11]
  
  out <- cbind(snpResult, fpcOut, inGao, inCurrantOuter, inCurrantInner) %>%
    as.data.table
  
  if(nrow(fpcOut) > 1){ print(paste(id,"fpc"))}
  if(nrow(inGao) > 1){ print(paste(id,"gao"))}
  if(nrow(inCurrantInner) > 1){ print(paste(id,"currantInner"))}
  if(nrow(inCurrantOuter) > 1){ print(paste(id,"currantOuter"))}
  
  return(out)
  
}) %>%
  rbindlist %>%
  .[, novel := ifelse(gaoSig=="N" & currantInnerSig=="N" & currantOuterSig=="N", "Y", "N")] %>% 
  .[vep, on = "ID"] %>%
  .[freqs, on = "ID"] %>%
  .[, EffAlleleFreq := ifelse(A1==ALT, ALT_FREQS, (1-ALT_FREQS))] %>%
  .[, EffAllele := A1] %>%
  .[, nonEffAllele := ifelse(A1==ALT, REF, ALT)] %>%
  .[, .(CHR, POS, ID, EffAllele, nonEffAllele, EffAlleleFreq, BETA, SE, P, pixel,
        nSNPsLocus, nPixelsLocus, fpcSig, topFPC, allFPCs, sameSentinel, FPCsentinel,
        FPCtopP, gaoSig, sameSentinelGao, currantOuterSig, sameSentinelCurrantOuter,
        currantInnerSig, sameSentinelCurrantInner, novel,  Consequence, SYMBOL,
        Existing_variation, DISTANCE, AF, PUBMED, PHENOTYPES)] %>%
  setnames(., "AF", "AF1000G")



## FPC results

vepFPC <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/misc/fpcSentinelsVEP20230328.txt") %>%
  .[, c("#Uploaded_variation", "Consequence", "SYMBOL", "Existing_variation", "DISTANCE", "AF", "PUBMED", "PHENOTYPES")] %>%
  setnames(., "#Uploaded_variation", "ID")

snps <- fpcSentinels[, ID]

fpcsAnnotated <- lapply(snps, function(id) {
  
  result <- fpcSentinels[ID==id]
  
  allLocusSNPs <-  c(result[,ID], result[,SNPsInLocus]) %>% 
    strsplit(., ",") %>%
    unlist
  
  pixelMatch <- sentinels[ID %in% allLocusSNPs]
  
  if(nrow(pixelMatch) > 0) {
    
    pixelOut <- data.table(pixelwiseSig = "Y",
                           toppixel = paste(pixelMatch[,pixel], collapse = ";"),
                           nPixels = paste(pixelMatch[, nPixelsLocus], collapse = ";"),
                           sameSentinel = ifelse(id %in% pixelMatch[,ID], "Y", "N"),
                           pixelSentinel = paste(pixelMatch[,ID], collapse = ";"),
                           pixelTopP = paste(pixelMatch[,P], collapse = ";"))
    
  } else {
    pixelOut <- data.table(pixelwiseSig = "N",
                           toppixel = NA,
                           nPixels = NA,
                           sameSentinel = NA,
                           pixelSentinel = NA,
                           pixelTopP = NA)
    
  }
  
  gaoMatch <- gao[rsid %in% allLocusSNPs]
  
  if(nrow(gaoMatch) > 0) {
    
    inGao <- data.table(gaoSig = "Y",
                        sameSentinelGao = ifelse(id %in% gaoMatch[,rsid], "Y", "N"))
    
  } else {
    inGao <- data.table(gaoSig = "N",
                        sameSentinelGao = NA)
  }
  
  currantInnerMatch <- currant1[SNP %in% allLocusSNPs]
  
  if(nrow(currantInnerMatch) > 0) {
    
    inCurrantInner <- data.table(currantInnerSig = "Y",
                                 sameSentinelCurrantInner = ifelse(id %in% currantInnerMatch[,SNP], "Y", "N"))
    
  } else {
    inCurrantInner <- data.table(currantInnerSig = "N",
                                 sameSentinelCurrantInner = NA)
  }
  
  currantOuterMatch <- currant2[SNP %in% allLocusSNPs]
  
  if(nrow(currantOuterMatch) > 0) {
    
    inCurrantOuter <- data.table(currantOuterSig = "Y",
                                 sameSentinelCurrantOuter = ifelse(any(currantOuterMatch[,SNP] %in% id), "Y", "N"))
    
  } else {
    inCurrantOuter <- data.table(currantOuterSig = "N",
                                 sameSentinelCurrantOuter = NA)
  }
  
  snpResult <- fpcSentinels[ID==id, -20]
  
  out <- cbind(snpResult, pixelOut, inGao, inCurrantOuter, inCurrantInner) %>%
    as.data.table
  
  if(nrow(pixelOut) > 1){ print(paste(id,"pixels"))}
  if(nrow(inGao) > 1){ print(paste(id,"gao"))}
  if(nrow(inCurrantInner) > 1){ print(paste(id,"currantInner"))}
  if(nrow(inCurrantOuter) > 1){ print(paste(id,"currantOuter"))}
  
  return(out)
  
}) %>%
  rbindlist %>%
  .[, novel := ifelse(gaoSig=="N" & currantInnerSig=="N" & currantOuterSig=="N", "Y", "N")] %>% 
  .[vepFPC, on = "ID"] %>%
  .[, EffAlleleFreq := A1_FREQ] %>%
  .[, EffAllele := A1] %>%
  .[, nonEffAllele := ifelse(A1==ALT, REF, ALT)] %>%
  .[, .(CHR, POS, ID, EffAllele, nonEffAllele, A1_FREQ, BETA, SE, P, FPC,
        nSNPsLocus, nFPCsLocus, clumpFPCs, pixelwiseSig, toppixel, nPixels, sameSentinel, pixelSentinel,
        pixelTopP, gaoSig, sameSentinelGao, currantOuterSig, sameSentinelCurrantOuter,
        currantInnerSig, sameSentinelCurrantInner, novel,  Consequence, SYMBOL,
        Existing_variation, DISTANCE, AF, PUBMED, PHENOTYPES)] %>%
  setnames(., "AF", "AF1000G") 



fwrite(pixelsAnnotated, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/annotatedSentinelsPixelwise.txt", sep = "\t")
fwrite(pixelsAnnotated[novel=="Y"], file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/annotatedSentinelsPixelwiseNovelOnly.txt", sep = "\t")

fwrite(fpcsAnnotated, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/sentinels/annotatedSentinelsFPCs.txt", sep = "\t")
fwrite(fpcsAnnotated[novel=="Y"], file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/sentinels/annotatedSentinelsFPCsNovelOnly.txt", sep = "\t")






pixelsAnnotated<- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/annotatedSentinelsPixelwise.txt")
fpcsAnnotated <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/sentinels/annotatedSentinelsFPCs.txt")

both <-  pixelsAnnotated[ID %in% fpcSentinelsAllSNPsBon, ID] %>% unique

upsetDataPixel <- pixelsAnnotated[P < (5E-8/29041), .(ID, fpcSig, gaoSig, currantInnerSig, currantOuterSig)] %>%
  .[, pixelWise := 1] %>%
  .[, FPC := ifelse(ID %in% fpcSentinelsAllSNPsBon, 1, 0)] %>%
  .[, Gao := ifelse(gaoSig=="Y", 1, 0)] %>%
  .[, CurrantInnerLayers := ifelse(currantInnerSig=="Y", 1, 0)] %>%
  .[, CurrantOuterLayers := ifelse(currantOuterSig=="Y", 1, 0)] %>%
  .[, .(ID, pixelWise, FPC, Gao, CurrantInnerLayers, CurrantOuterLayers)]

upsetFPCOnly <- fpcsAnnotated[!(ID %in% sentinelsAllSNPsBon) & P < (5E-8/6), .(ID, pixelwiseSig, gaoSig, currantInnerSig, currantOuterSig)] %>%
  .[, pixelWise := 0] %>%
  .[, FPC := 1] %>%
  .[, Gao := ifelse(gaoSig=="Y", 1, 0)] %>%
  .[, CurrantInnerLayers := ifelse(currantInnerSig=="Y", 1, 0)] %>%
  .[, CurrantOuterLayers := ifelse(currantOuterSig=="Y", 1, 0)] %>%
  .[, .(ID, pixelWise, FPC, Gao, CurrantInnerLayers, CurrantOuterLayers)]

upsetData <- rbind(upsetDataPixel, upsetFPCOnly)


png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/plots/upsetAllLoci.png", width = 1200, height = 600)
upset(upsetData) %>% print
dev.off()























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

