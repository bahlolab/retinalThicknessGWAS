library(data.table)
library(magrittr)
library(tidyverse)
library(patchwork)
library(UpSetR)
library(stringr)

## exploring overlap of fPC and pixel-wise approaches
## and comparing loci with previous hits

fpcSentinels <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/sentinels/allChr_sentinel_clumpThresh0.001_withOverlap.txt") %>%
  .[nSNPsLocus >= 5]

sentinels <- lapply(c(1:22, "X"), function(chr) {
  
  chrSent <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/chr",chr,"sentinels_clumpThresh0.001_withOverlap.txt") %>%
    fread() %>%
    .[nSNPsLocus >= 5]
  return(chrSent)
}) %>%
  rbindlist(., idcol = "CHR")


## upsettr plots showing overlap
fpcSentinelsAllSNPs <- c(fpcSentinels[nSNPsLocus >= 5,ID], fpcSentinels[nSNPsLocus >= 5,SNPsInLocus]) %>% 
  strsplit(., ",") %>%
  unlist

fpcSentinelsAllSNPsBon <- c(fpcSentinels[nSNPsLocus >= 5 & P < (5E-8/6),ID], fpcSentinels[nSNPsLocus >= 5 & P < (5E-8/6),SNPsInLocus]) %>% 
  strsplit(., ",") %>%
  unlist

sentinelsAllSNPs <- c(sentinels[nSNPsLocus >= 5,ID], sentinels[nSNPsLocus >= 5,SNPsInLocus]) %>% 
  strsplit(., ",") %>%
  unlist

sentinelsAllSNPsBon <- c(sentinels[nSNPsLocus >= 5 & P < (5E-8/29041),ID], sentinels[nSNPsLocus >= 5 & P < (5E-8/29041),SNPsInLocus]) %>% 
  strsplit(., ",") %>%
  unlist

pixelsAnnotated<- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/annotatedSentinelsPixelwise.txt")
fpcsAnnotated <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/annotatedSentinelsFPCs.txt")

both <-  pixelsAnnotated[nSNPsLocus >= 5 & P < (5E-8/29041) & ID %in% fpcSentinelsAllSNPsBon, ID] %>% unique

upsetDataPixel <- pixelsAnnotated[nSNPsLocus >= 5 & P < (5E-8/29041), .(ID, fpcSig, gaoSig, currantInnerSig, currantOuterSig, zekavatSig)] %>%
  unique %>%
  .[, pixelWise := 1] %>%
  .[, FPC := ifelse(ID %in% fpcSentinelsAllSNPsBon, 1, 0)] %>%
  .[, Gao := ifelse(gaoSig=="Y", 1, 0)] %>%
  .[, CurrantInnerLayers := ifelse(currantInnerSig=="Y", 1, 0)] %>%
  .[, CurrantOuterLayers := ifelse(currantOuterSig=="Y", 1, 0)] %>%
  .[, Zekavat := ifelse(zekavatSig=="Y", 1, 0)] %>%
  .[, .(ID, pixelWise, FPC, Gao, CurrantInnerLayers, CurrantOuterLayers, Zekavat)]

upsetFPCOnly <- fpcsAnnotated[!(ID %in% sentinelsAllSNPsBon) & nSNPsLocus >= 5  & P < (5E-8/6), .(ID, pixelwiseSig, gaoSig, currantInnerSig, currantOuterSig, zekavatSig)] %>%
  unique %>%
  .[, pixelWise := 0] %>%
  .[, FPC := 1] %>%
  .[, Gao := ifelse(gaoSig=="Y", 1, 0)] %>%
  .[, CurrantInnerLayers := ifelse(currantInnerSig=="Y", 1, 0)] %>%
  .[, CurrantOuterLayers := ifelse(currantOuterSig=="Y", 1, 0)] %>%
  .[, Zekavat := ifelse(zekavatSig=="Y", 1, 0)] %>%
  .[, .(ID, pixelWise, FPC, Gao, CurrantInnerLayers, CurrantOuterLayers, Zekavat)]

upsetData <- rbind(upsetDataPixel, upsetFPCOnly)


png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/plots/upsetAllLoci.png", width = 1200, height = 600)
upset(upsetData, nsets = 6, text.scale = 2) %>% print
dev.off()




## How many of the known loci aren't identified through our analyses?

## read in previous SNPs
rsids <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/misc/rsids_Gao.txt") %>%
  .[, start := NULL]
gao <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/misc/macularThicknessLoci_Gao.txt") %>%
  rsids[., on = c("chr", "pos")]

currant1 <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/misc/retinalThicknessLoci_Currant.csv")
currant2 <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/misc/retinalThicknessLoci_Currant_2023.txt")

zekavat <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/misc/Zekavat/genomeWideSigZekavat2024.csv") 

gao[,rsid] %in% c(sentinelsAllSNPsBon, fpcSentinelsAllSNPsBon) %>% table
currant1[,SNP] %in% c(sentinelsAllSNPsBon, fpcSentinelsAllSNPsBon) %>% table
currant2[,SNP] %in% c(sentinelsAllSNPsBon, fpcSentinelsAllSNPsBon) %>% table
unique(zekavat[,rsid]) %in% c(sentinelsAllSNPsBon, fpcSentinelsAllSNPsBon) %>% table



gaoCompare <- data.table(SNP = gao[,rsid],
                         MAF = ifelse(gao[,AF] < 0.5, gao[,AF], 1-gao[,AF]),
                         log10P = (-1)*log(gao[,p], 10),
                         RTsigBon = ifelse(gao[,rsid] %in% c(sentinelsAllSNPsBon, fpcSentinelsAllSNPsBon), "Y", "N"),
                         RTsigGW = ifelse(gao[,rsid] %in% c(sentinelsAllSNPs, fpcSentinelsAllSNPs), "Y", "N"))

ggplot(gaoCompare, aes(x=MAF, y=log10P , colour=RTsigBon)) +
  geom_point()

ggplot(gaoCompare, aes(x=MAF, y=log10P , colour=RTsigGW)) +
  geom_point()




currant1compare <- data.table(SNP = currant1[,SNP],
                              MAF = ifelse(currant1[,AF] < 0.5, currant1[,AF], 1-currant1[,AF]),
                              MTAGlog10P = (-1)*log(unlist(currant1[,c("MTAG p-value")]), 10),
                              RNFLlog10P = (-1)*log(unlist(currant1[,"RNFL p-value"]), 10),
                              GCIPLlog10P = (-1)*log(unlist(currant1[,"GCIPL p-value"]), 10),
                              RTsigBon = ifelse(currant1[,SNP] %in% c(sentinelsAllSNPsBon, fpcSentinelsAllSNPsBon), "Y", "N"),
                              RTsigGW = ifelse(currant1[,SNP] %in% c(sentinelsAllSNPs, fpcSentinelsAllSNPs), "Y", "N"))

currant1compareLongCommon <- na.omit(currant1compare[MAF > 0.05]) %>%
  melt(., id.vars = c("SNP", "MAF", "RTsig"), 
       measure.vars = c("MTAGlog10P", "RNFLlog10P", "GCIPLlog10P"), 
       variable.name = "analysis", 
       value.name = "log10P")

ggplot(currant1compare, aes(x=MAF, y=MTAGlog10P, colour=RTsig)) +
  geom_point()

ggplot(currant1compareLongCommon, aes(x=SNP, y=log10P, colour=RTsig, shape = analysis)) +
  geom_point()


allKnownSNPs <- data.table(ID = c(gao[,rsid], currant1[,SNP], currant2[,SNP], unique(zekavat[,rsid])))
fwrite(allKnownSNPs, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/misc/allKnownSNPs.txt")




## previously identified loci - pixel wise results - minP comparison.
pixelResultsKnownLoci <- lapply(c(1:22), function(chr) {
  
  dir <- paste0("/vast/scratch/users/jackson.v/retThickness/GWASfollowup/pixelWiseResultsKnownLoci/results/chr",chr)
  print(paste("chr",chr))
  
  pixResults <- lapply(c(1:119), function(slice) {
    
    
    file <- paste0(dir,"/chr",chr,"Slice",slice,"_sentinels.txt")
    
    if(file.exists(file)){
      sliceResults <- fread(file) %>%
        setnames(., "#POS", "POS") %>%
        .[, chrom := chr] %>%
        .[, .(chrom, POS, ID, BETA, P, pixel)] %>%
        setnames(., c("chr", "pos", "ID", "beta", "P", "Trait"))
      
      return(sliceResults)
    }
    
  }) %>%
    rbindlist 
  
  fpcRes <- lapply(c(1:6), function(i) {
    print(paste(i))
    
    paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/results/chr",chr,"/chr",chr,"EUR.fpc",i,".glm.linear") %>%
      fread(., select = c("#CHROM", "POS", "ID", "BETA", "P")) %>%
      .[, Trait := paste0("FPC",i)] %>%
      .[ID %in% pixResults[,ID]] %>%
      setnames(., c("chr", "pos", "ID", "beta", "P", "Trait")) %>%
      return(.)
  }) %>%
    rbindlist 
  
  
  chrResults <- rbind(pixResults, fpcRes) %>%
    .[, .SD[which.min(as.numeric(P))], by = ID] 
  
  return(chrResults)
  
}) %>%
  rbindlist



freqs <- lapply(c(1:22), function(chr) {
  print(paste("chr",chr))
  
  paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/results/chr",chr,"/chr",chr,"EUR.fpc1.glm.linear") %>%
    fread(., select = c("ID", "A1_FREQ")) %>%
    .[ID %in% pixelResultsKnownLoci[,ID]] %>%
    .[, A1_FREQ := ifelse(A1_FREQ < 0.5, A1_FREQ, (1 - A1_FREQ))] %>%
    setnames(., c("ID", "MAF")) %>%
    return(.)
}) %>%
  rbindlist 

knownLoci <- pixelResultsKnownLoci[,pos := as.numeric(pos)] %>%
  .[freqs, on = c("ID")]

gaoCompare <- gao[knownLoci, on  = c("chr", "pos", "ID")] %>%
  .[, log10pGao := (-1)*log(p, 10)] %>%
  .[, log10pRT := (-1)*log(as.numeric(P), 10)] 


png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/gaoCompareAll.png", width = 500, height = 500)
ggplot(gaoCompare, aes(x=log10pGao, y=log10pRT, size = MAF)) + 
  geom_point()  +
  geom_abline(intercept=0, slope = 1, linetype = "dashed") +
  theme_bw() 
dev.off()

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/gaoCompareNotSig.png", width = 500, height = 500)
ggplot(gaoCompare[as.numeric(P) >= 5E-8/29041], aes(x=log10pGao, y=log10pRT, size = MAF)) + 
  geom_point()  +
  geom_abline(intercept=0, slope = 1, linetype = "dashed") +
  geom_hline(yintercept = (-1)*log10(5E-8), linetype = "dotted", color = "red") +
  theme_bw() 
dev.off()

currant1Compare <- currant1[knownLoci, on  = c("chr", "pos", "ID")] %>%
  .[, log10pCurrant := (-1)*log(p, 10)] %>%
  .[, log10pRT := (-1)*log(as.numeric(P), 10)] 

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/currant1CompareAll.png", width = 500, height = 500)
ggplot(currant1Compare, aes(x=log10pCurrant, y=log10pRT, size = MAF)) + 
  geom_point()  +
  geom_abline(intercept=0, slope = 1, linetype = "dashed") +
  theme_bw() 
dev.off()

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/currant1CompareNotSig.png", width = 500, height = 500)
ggplot(currant1Compare[as.numeric(P) >= 5E-8/29041], aes(x=log10pCurrant, y=log10pRT, size = MAF)) + 
  geom_point()  +
  geom_abline(intercept=0, slope = 1, linetype = "dashed") +
  geom_hline(yintercept = (-1)*log10(5E-8), linetype = "dotted", color = "red") +
  theme_bw() 
dev.off()


currant2Compare <- currant2[knownLoci, on  = c("chr", "pos", "ID")] %>%
  .[, log10pCurrant := (-1)*log(p, 10)] %>%
  .[, log10pRT := (-1)*log(as.numeric(P), 10)] 

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/currant2CompareAll.png", width = 500, height = 500)
ggplot(currant2Compare, aes(x=log10pCurrant, y=log10pRT, size = MAF)) + 
  geom_point()  +
  geom_abline(intercept=0, slope = 1, linetype = "dashed") +
  theme_bw() 
dev.off()

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/currant2CompareNotSig.png", width = 500, height = 500)
ggplot(currant2Compare[as.numeric(P) >= 5E-8/29041], aes(x=log10pCurrant, y=log10pRT, size = MAF)) + 
  geom_point()  +
  geom_abline(intercept=0, slope = 1, linetype = "dashed") +
  geom_hline(yintercept = (-1)*log10(5E-8), linetype = "dotted", color = "red") +
  theme_bw() 
dev.off()


zekavatCompare <- zekavat[knownLoci, on  = c("rsid" ="ID")] %>%
  .[, .SD[which.min(as.numeric(p_value.x))], by = rsid] %>%
  .[, log10pZekavat := (-1)*log(p_value.x, 10)] %>%
  .[, log10pRT := (-1)*log(as.numeric(P), 10)] 

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/zekavatCompareAll.png", width = 500, height = 500)
ggplot(zekavatCompare, aes(x=log10pZekavat, y=log10pRT, size = MAF)) + 
  geom_point()  +
  geom_abline(intercept=0, slope = 1, linetype = "dashed") +
  theme_bw() 
dev.off()

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/zekavatCompareNotSig.png", width = 500, height = 500)
ggplot(zekavatCompare[as.numeric(P) >= 5E-8/29041], aes(x=log10pZekavat, y=log10pRT, size = MAF)) + 
  geom_point()  +
  geom_abline(intercept=0, slope = 1, linetype = "dashed") +
  geom_hline(yintercept = (-1)*log10(5E-8), linetype = "dotted", color = "red") +
  theme_bw() 
dev.off()













## read in pixel-wise reasults for all FPC sentinels
## previously identified loci - pixel wise results - minP comparison.
pixelResultsFPCloci <- lapply(c(1:22), function(chr) {
  
  dir <- paste0("/vast/scratch/users/jackson.v/retThickness/GWASfollowup/pixelWiseResultsFPCsentinels/results/chr",chr)
  print(paste("chr",chr))
  
  pixResults <- lapply(c(1:119), function(slice) {
    
    
    file <- paste0(dir,"/chr",chr,"Slice",slice,"_sentinels.txt")
    
    sliceResults <- fread(file) %>%
      setnames(., "#POS", "POS") %>%
      .[, chrom := chr] %>%
      .[, .(chrom, POS, ID, BETA, P, pixel)] %>%
      setnames(., c("chr", "pos", "ID", "beta", "P", "Trait"))
    
    return(sliceResults)
    
  }) %>%
    rbindlist 
  
  
  return(pixResults)
  
}) %>%
  rbindlist

# find min P for each ID
minPixP <- pixelResultsFPCloci %>%
  .[, .SD[which.min(as.numeric(P))], by = ID] 


minPixP[order(P)] %>% head























## Now examine overlap/differences between fPC and pixel-wise approaches..
## average patterns - those sig in fPC analysis; and those not

fpcSig <- results[ID %in% sentinelsAllSNPsBon & ID %in% fpcSentinelsAllSNPsBon] %>%
  .[, mean(.SD, na.rm = T), by = pixel, .SDcols = c("BETA", "P")]

fpcSigBetas <- results[ID %in% sentinelsAllSNPsBon & ID %in% fpcSentinelsAllSNPsBon] %>%
  .[, mean(BETA, na.rm = T), by = pixel] %>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert = T)] %>%
  .[!is.na(y)] 

fpcSigP <- results[ID %in% sentinelsAllSNPsBon & ID %in% fpcSentinelsAllSNPsBon] %>%
  .[, mean(P, na.rm = T), by = pixel] %>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert = T)] %>%
  .[!is.na(y)] %>%
  .[, log10P := (-1)*log(V1, 10)]


fpcSigNP <- results[ID %in% sentinelsAllSNPsBon & ID %in% fpcSentinelsAllSNPsBon] %>%
  .[P < 5E-8/29041, .N, by = pixel] %>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert = T)] %>%
  .[!is.na(y)] 


ggplot(fpcSigBetas) +
  geom_tile(aes_string(x = "x", y = "y", fill = "V1")) +
  # geom_path(aes(x=col,y=row,color=area),data = areas, size = 0.5) +
  scale_fill_gradient2() +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "bottom") +
  ggtitle("Mean effect - SNPs significant in fPC analysis")

ggplot(fpcSigP) +
  geom_tile(aes_string(x = "x", y = "y", fill = "log10P")) +
  # geom_path(aes(x=col,y=row,color=area),data = areas, size = 0.5) +
  scale_fill_gradient2() +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "bottom") +
  ggtitle("Mean -log10P - SNPs significant in fPC analysis")

ggplot(fpcSigNP) +
  geom_tile(aes_string(x = "x", y = "y", fill = "N")) +
  # geom_path(aes(x=col,y=row,color=area),data = areas, size = 0.5) +
  scale_fill_gradient2() +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "bottom") +
  ggtitle("Count of significant associations at pixel - SNPs significant in fPC analysis")



fpcNotSigBetas <- results[ID %in% sentinelsAllSNPsBon & !(ID %in% fpcSentinelsAllSNPsBon)] %>%
  .[, mean(BETA, na.rm = T), by = pixel] %>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert = T)] %>%
  .[!is.na(y)] 

fpcNotSigP <- results[ID %in% sentinelsAllSNPsBon & !(ID %in% fpcSentinelsAllSNPsBon)] %>%
  .[, mean(P, na.rm = T), by = pixel] %>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert = T)] %>%
  .[!is.na(y)] %>%
  .[, log10P := (-1)*log(V1, 10)]

fpcNotSigNP <- results[ID %in% sentinelsAllSNPsBon & !(ID %in% fpcSentinelsAllSNPsBon)] %>%
  .[P < 5E-8/29041, .N, by = pixel] %>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert = T)] %>%
  .[!is.na(y)] 

ggplot(fpcNotSigBetas) +
  geom_tile(aes_string(x = "x", y = "y", fill = "V1")) +
  # geom_path(aes(x=col,y=row,color=area),data = areas, size = 0.5) +
  scale_fill_gradient2() +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "bottom") +
  ggtitle("Mean effect - SNPs not significant in fPC analysis")

ggplot(fpcNotSigP) +
  geom_tile(aes_string(x = "x", y = "y", fill = "log10P")) +
  # geom_path(aes(x=col,y=row,color=area),data = areas, size = 0.5) +
  scale_fill_gradient2() +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "bottom") +
  ggtitle("Mean -log10P - SNPs not significant in fPC analysis")

ggplot(fpcNotSigNP) +
  geom_tile(aes_string(x = "x", y = "y", fill = "N")) +
  # geom_path(aes(x=col,y=row,color=area),data = areas, size = 0.5) +
  scale_fill_gradient2() +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "bottom") +
  ggtitle("Count of significant associations at pixel - SNPs not significant in fPC analysis")







## average patterns - those sig in fPC analysis; and those not
pixelSig <- fpcResults[ID %in% sentinelsAllSNPsBon & ID %in% fpcSentinelsAllSNPsBon] %>%
  melt(., id.vars = "ID", measure.vars = paste0("P_FPC", c(1:6)) , variable.name = "FPC", value.name = "P") %>%
  unique %>%
  .[P < 5E-8/6, .N, by = FPC] 

pixelNotSig <- fpcResults[!(ID %in% sentinelsAllSNPsBon) & ID %in% fpcSentinelsAllSNPsBon] %>%
  melt(., id.vars = "ID", measure.vars = paste0("P_FPC", c(1:6)) , variable.name = "FPC", value.name = "P") %>%
  unique %>%
  .[P < 5E-8/6, .N, by = FPC] 



## Examine overlap/differences between fPC and pixel-wise approaches

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





