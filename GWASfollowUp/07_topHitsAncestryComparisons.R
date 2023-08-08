library(data.table)
library(magrittr)
library(tidyverse)
library(patchwork)
library(scales)

pixels <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/pixels.txt") %>%
setnames(., c("pixel", "y", "x"))

sents <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/annotatedSentinelsPixelwise.txt") %>%
    .[, CHROM := ifelse(CHR==23, "X", CHR)] %>%
  .[, SE := abs(BETA / qnorm(P/2))] %>%
  .[, .(CHROM, POS, ID, EffAllele, EffAlleleFreq, BETA, SE, P, BonferroniSig, topPixel)] %>%
  setnames(., c("CHROM", "POS", "ID", "A1", "A1_FREQ", "BETA", "SE", "P","BonferroniSig", "pixel" ))

sentsAFR <- lapply(c(1:22, "X"), function(chr) {

    print(paste(chr))


    ANC <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinelsAFR/chr",chr,"EURsentinels_AFRassocsAllPixels.txt") %>%
    fread()

    chrResult <- lapply(sents[CHROM==as.character(chr),ID], function(snp) {

        topPix <- sents[ID==snp, pixel]

        resultANC <- ANC[ID==snp & pixel==topPix]

        return(resultANC)

    }) %>%
    rbindlist %>%
    .[, c("CHROM", "POS", "ID", "A1", "A1_FREQ", "BETA", "SE", "P", "pixel" )]

    return(chrResult)}) %>%

rbindlist %>%
.[A1_FREQ < 0.001, c("BETA", "P") := NA]


sentsCSA <- lapply(c(1:22, "X"), function(chr) {

    print(paste(chr))

    ANC <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinelsCSA/chr",chr,"EURsentinels_CSAassocsAllPixels.txt") %>%
    fread()

    chrResult <- lapply(sents[,ID], function(snp) {

        topPix <- sents[ID==snp, pixel]

        resultANC <- ANC[ID==snp & pixel==topPix]

        return(resultANC)

    }) %>%

    rbindlist %>% 
    .[, c("CHROM", "POS", "ID", "A1", "A1_FREQ", "BETA", "SE", "P", "pixel" )]

    return(chrResult)}) %>%

rbindlist %>%
.[A1_FREQ < 0.001, c("BETA", "P") := NA]

allSents <- sents %>%
    merge(., sentsAFR, by = c("CHROM", "POS", "ID", "pixel"), suffixes = c("", ".AFR"), all = T) %>%
    merge(., sentsCSA, on = c("POS", "ID", "REF", "ALT", "pixel"), suffixes = c("", ".CSA"), all =T) %>%
  setnames(., old = c("A1_FREQ", "BETA", "SE", "P"), new = c("A1_FREQ.EUR", "BETA.EUR", "SE.EUR", "P.EUR")) %>%
  .[, BETA.AFR := case_when(A1==A1.AFR ~ BETA.AFR,
                            T ~ (-1)*BETA.AFR)] %>%
   .[, BETA.CSA := case_when(A1==A1.CSA ~ BETA.CSA,
                            T ~ (-1)*BETA.CSA)] %>%
   .[, A1_FREQ.AFR := case_when(A1==A1.AFR ~ A1_FREQ.AFR,
                            T ~ (-1)*A1_FREQ.AFR)] %>%
   .[, A1_FREQ.CSA := case_when(A1==A1.CSA ~ A1_FREQ.CSA,
                            T ~ (-1)*A1_FREQ.CSA)] %>%
    .[, EurAfrConcordant := case_when(sign(BETA.EUR)==sign(BETA.AFR) ~ 1, 
                                      sign(BETA.EUR)!=sign(BETA.AFR) ~ 0,
                                      is.na(BETA.AFR) ~ NA_real_ )]  %>%
    .[, EurCsaConcordant := case_when(sign(BETA.EUR)==sign(BETA.CSA) ~ 1, 
                                      sign(BETA.EUR)!=sign(BETA.CSA) ~ 0,
                                      is.na(BETA.CSA) ~ NA_real_ )]                



allSentsBonf <- allSents[BonferroniSig == "Y"] %>%
   .[, .(CHROM, POS, ID, pixel, A1, A1_FREQ.EUR, BETA.EUR, SE.EUR, P.EUR, A1_FREQ.AFR, BETA.AFR, SE.AFR, P.AFR, A1_FREQ.CSA, BETA.CSA, SE.CSA, P.CSA, EurAfrConcordant, EurCsaConcordant)]          

thresh <- 0.05/nrow(allSentsBonf)
threshLabel <- scientific(thresh, digits = 3)

# Create scatter plot of BETA vs BETA.AFR, colored by P.AFR threshold
afrPlot <- ggplot(allSentsBonf, aes(x = BETA.EUR, y = BETA.AFR, color = P.AFR < thresh)) +
  geom_point() + 
  geom_errorbar(aes(ymin = BETA.AFR - SE.AFR , ymax = BETA.AFR + SE.AFR)) + 
  geom_errorbarh(aes(xmin = BETA.EUR - SE.EUR, xmax = BETA.EUR + SE.EUR))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
  scale_color_manual(values = c("TRUE" = "#CC0066", "FALSE" = "#6699CC")) +
  labs(x = "BETA (EUR)", y = "BETA (AFR)", color =  substitute(paste("P"[AFR], "<", threshLabel), list(threshLabel = threshLabel)))

# Create scatter plot of BETA vs BETA.CSA, colored by P.CSA threshold
csaPlot <- ggplot(allSentsBonf, aes(x = BETA.EUR, y = BETA.CSA, color = P.CSA < thresh)) +
  geom_point() +
  geom_errorbar(aes(ymin = BETA.CSA - SE.CSA , ymax = BETA.CSA + SE.CSA)) + 
  geom_errorbarh(aes(xmin = BETA.EUR - SE.EUR, xmax = BETA.EUR + SE.EUR))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
  scale_color_manual(values = c("TRUE" = "#CC0066", "FALSE" = "#6699CC")) +
  labs(x = "BETA (EUR)", y = "BETA (CSA)", color =  substitute(paste("P"[CSA], "<", threshLabel), list(threshLabel = threshLabel)))

  png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/plots/ancestriesComparisonBonfSigSentinels.png", width = 400, height = 800)
  afrPlot / csaPlot %>%
  print
  dev.off()


  
 ## Look at correlations of scan-wise effects for each SNP.
  results <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/allSentinelsAllPixelsResults_clumpThresh0.001_withOverlap.csv")
  

  # corrEURAFR <- lapply(c(1:22, "X"), function(chr) {
    corrEURAFR <- lapply(c(1:22), function(chr) {
      
    print(paste(chr))
    
    
    ANC <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinelsAFR/chr",chr,"EURsentinels_AFRassocsAllPixels.txt") %>%
      fread()
    
    chrResult <- lapply(sents[CHROM==as.character(chr),ID], function(snp) {
      
      ID <- snp
      A1_FREQ <- ANC[ID==snp, A1_FREQ] %>% unique()
      
      compDT <- results[ID==snp, .(pixel, A1, BETA)] %>%
        .[ANC[ID==snp, .(pixel, A1, BETA)], on = "pixel"]
      
      if(compDT[1,A1] == compDT[1,i.A1]) {
        
        corr <- cor(compDT[,BETA], compDT[,i.BETA], method = "spearman")
        
      } else {
        print(paste("alleles switched -",ID))
        corr <- cor(compDT[,BETA], (-1)*compDT[,i.BETA], method = "spearman")
      }
      
      resultANC <- data.table(ID, A1_FREQ, corr)
      
      return(resultANC)
      
    }) %>%
      rbindlist 
    
    return(chrResult)}) %>%
    
    rbindlist %>%
    .[A1_FREQ >= 0.001]
  
  

  
    corrEURCSA <- lapply(c(1:22), function(chr) {
      # corrEURCSA <- lapply(c(1:22, "X"), function(chr) {
    
    print(paste(chr))
    
    
    ANC <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinelsCSA/chr",chr,"EURsentinels_CSAassocsAllPixels.txt") %>%
      fread()
    
    chrResult <- lapply(sents[CHROM==as.character(chr),ID], function(snp) {
      
      ID <- snp
      A1_FREQ <- ANC[ID==snp, A1_FREQ] %>% unique()
      
      compDT <- results[ID==snp, .(pixel, A1, BETA)] %>%
        .[ANC[ID==snp, .(pixel, A1, BETA)], on = "pixel"]
      
      if(compDT[1,A1] == compDT[1,i.A1]) {
        
        corr <- cor(compDT[,BETA], compDT[,i.BETA], method = "spearman")
        
      } else {
        print(paste("alleles switched -",ID))
      corr <- cor(compDT[,BETA], (-1)*compDT[,i.BETA], method = "spearman")
      }
      
      resultANC <- data.table(ID, A1_FREQ, corr)
      
      return(resultANC)
      
    }) %>%
      rbindlist 
    
    return(chrResult)}) %>%
    
    rbindlist %>%
    .[A1_FREQ >= 0.001]
  
  
 corrAFR <- ggplot(corrEURAFR, aes(x = corr)) +
   geom_density() +
   geom_vline(xintercept = mean(corrEURAFR[,corr]), linetype = "dashed", color = "grey") +
   labs(x = "Spearman's rho") +
   theme_minimal() +
   xlim(-1, 1)
 
 corrCSA <- ggplot(corrEURCSA, aes(x = corr)) +
   geom_density() +
   geom_vline(xintercept = corrEURCSA[,corr], linetype = "dashed", color = "grey") +
   labs(x = "Spearman's rho") +
   theme_minimal() +
   xlim(-1, 1)

 
 png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/plots/ancestriesComparisonScanwiseCorrs.png", height = 800, width = 400)
 corrAFR / corrCSA %>%
   print
 dev.off()

 
 corrAFR <- ggplot(corrEURAFR[ID %in% sents[BonferroniSig == "Y", ID]], aes(x = corr)) +
   geom_density(color = "#CC0066") +
   geom_vline(xintercept = median(corrEURAFR[ID %in% sents[BonferroniSig == "Y",  ID], corr], na.rm=T), linetype = "dashed", color = "grey") +
   labs(x = "Correlation of pixel effects scan-wide") +
   theme_minimal() +
   xlim(-1, 1)
 
 corrCSA <- ggplot(corrEURCSA[ID %in% sents[BonferroniSig == "Y", ID] ], aes(x = corr)) +
   geom_density(color = "#CC0066") +
   geom_vline(xintercept = median(corrEURCSA[ID %in% sents[BonferroniSig == "Y", ID], corr], na.rm=T), linetype = "dashed", color = "grey") +
   labs(x = "Correlation of pixel effects scan-wide") +
   theme_minimal() +
   xlim(-1, 1)
 
 
 png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/plots/ancestriesComparisonScanwiseCorrsBonfSig.png", height = 800, width = 400)
 corrAFR / corrCSA %>%
   print
 dev.off()
 
 corrEURAFR[!is.na(corr)] %>% nrow
 corrEURAFR[corr > 0.5] %>% nrow
 
 corrEURAFR[ID %in% sents[BonferroniSig == "Y", ID] & !is.na(corr)] %>% nrow
 corrEURAFR[ID %in% sents[BonferroniSig == "Y", ID] &corr > 0.5] %>% nrow
 
 corrEURAFR[ID %in% sents[BonferroniSig == "Y",ID], corr] %>% median(., na.rm=T)

 corrEURCSA[!is.na(corr)] %>% nrow
 corrEURCSA[corr > 0.5] %>% nrow
 
 corrEURCSA[ID %in% sents[BonferroniSig == "Y", ID] & !is.na(corr)] %>% nrow
 corrEURCSA[ID %in% sents[BonferroniSig == "Y", ID] &corr > 0.5] %>% nrow

 corrEURCSA[ID %in% sents[BonferroniSig == "Y",ID], corr] %>% median(., na.rm=T)
 
# thresh <- 0.05/nrow(allSents)

# # Create scatter plot of BETA vs BETA.AFR, colored by P.AFR threshold
# afrPlot <- ggplot(allSents, aes(x = BETA, y = BETA.AFR, color = P.AFR < thresh)) +
#   geom_point() +
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
#   scale_color_manual(values = c("TRUE" = "#CC0066", "FALSE" = "#6699CC")) +
#   labs(x = "BETA (EUR)", y = "BETA (AFR)", color = bquote(paste("P[AFR] <", scientific(thresh, digits = 3))))

# # Create scatter plot of BETA vs BETA.CSA, colored by P.CSA threshold
# csaPlot <- ggplot(allSents, aes(x = BETA, y = BETA.CSA, color = P.CSA < thresh)) +
#   geom_point() +
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
#   scale_color_manual(values = c("TRUE" = "#CC0066", "FALSE" = "#6699CC")) +
#   labs(x = "BETA", y = "BETA.CSA", color = paste("P.CSA <", scientific(thresh, digits = 3)))

#   png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/plots/ancestriesComparisonGWSigSentinels.png", width = 1200, height = 600)
#   afrPlot + csaPlot %>%
#   print
#   dev.off()


fwrite(allSentsBonf, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/ancestryComparisons.csv")


allSentsBonf[, EurAfrConcordant] %>% table
allSentsBonf[, EurCsaConcordant] %>% table

cor.test(allSentsBonf[, BETA.EUR],  allSentsBonf[, BETA.AFR])
cor.test(allSentsBonf[, BETA.EUR],  allSentsBonf[, BETA.CSA])







## FPC results
sents <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/annotatedSentinelsFPCs.txt") %>%
    .[, CHROM := ifelse(CHR==23, "X", CHR)] %>%
    .[EffAlleleFreq > 0.001] %>%
    .[, SE := abs(BETA / qnorm(P/2))] %>%
    .[, .(CHROM, POS, ID, EffAllele, EffAlleleFreq, BETA, SE, P, BonferroniSig, topFPC)] %>%
    setnames(., c("CHROM", "POS", "ID", "A1", "A1_FREQ", "BETA", "SE", "P","BonferroniSig", "FPC" ))

    ANC <-  fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/AFRsentinels/allChr_EURsentinels_AFRassocsAllFPCs.txt")

    sentsAFR <- lapply(sents[,ID], function(snp) {

        topFPC <- sents[ID==snp, FPC]

        resultANC <- ANC[ID==snp & FPC==topFPC]

        return(resultANC)

    }) %>%
    rbindlist %>%
    .[, c("CHROM", "POS", "ID", "A1", "A1_FREQ", "BETA", "SE", "P", "FPC" )]%>%
  .[A1_FREQ < 0.001, c("BETA", "P") := NA]


    ANC <-  fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/CSAsentinels/allChr_EURsentinels_CSAassocsAllFPCs.txt")

    sentsCSA <- lapply(sents[,ID], function(snp) {

        topFPC <- sents[ID==snp, FPC]

        resultANC <- ANC[ID==snp & FPC==topFPC]

        return(resultANC)

    }) %>%
    rbindlist %>%
    .[, c("CHROM", "POS", "ID", "A1", "A1_FREQ", "BETA", "SE", "P", "FPC" )]%>%
  .[A1_FREQ < 0.001, c("BETA", "P") := NA]



allSents <- sents %>%
    merge(., sentsAFR, by = c("CHROM", "POS", "ID", "FPC"), suffixes = c("", ".AFR"), all = T) %>%
    merge(., sentsCSA, on = c("POS", "ID", "REF", "ALT", "FPC"), suffixes = c("", ".CSA"), all =T) %>%
    setnames(., old = c("A1_FREQ", "BETA", "SE", "P"), new = c("A1_FREQ.EUR", "BETA.EUR", "SE.EUR", "P.EUR")) %>%
   .[, BETA.AFR := case_when(A1==A1.AFR ~ BETA.AFR,
                            T ~ (-1)*BETA.AFR)] %>%
   .[, BETA.CSA := case_when(A1==A1.CSA ~ BETA.CSA,
                            T ~ (-1)*BETA.CSA)] %>%
   .[, A1_FREQ.AFR := case_when(A1==A1.AFR ~ A1_FREQ.AFR,
                            T ~ (-1)*A1_FREQ.AFR)] %>%
   .[, A1_FREQ.CSA := case_when(A1==A1.CSA ~ A1_FREQ.CSA,
                            T ~ (-1)*A1_FREQ.CSA)] %>%
    .[, EurAfrConcordant := case_when(sign(BETA.EUR)==sign(BETA.AFR) ~ 1, 
                                      sign(BETA.EUR)!=sign(BETA.AFR) ~ 0,
                                      is.na(BETA.AFR) ~ NA_real_ )]  %>%
    .[, EurCsaConcordant := case_when(sign(BETA.EUR)==sign(BETA.CSA) ~ 1, 
                                      sign(BETA.EUR)!=sign(BETA.CSA) ~ 0,
                                      is.na(BETA.CSA) ~ NA_real_ )]                



allSentsBonf <- allSents[BonferroniSig == "Y"] %>%
   .[, .(CHROM, POS, ID, FPC, A1, A1_FREQ.EUR, BETA.EUR, SE.EUR, P.EUR, A1_FREQ.AFR, BETA.AFR, SE.AFR, P.AFR, A1_FREQ.CSA, BETA.CSA, SE.CSA, P.CSA, EurAfrConcordant, EurCsaConcordant)]          

thresh <- 0.05/nrow(allSentsBonf)
threshLabel <- scientific(thresh, digits = 3)

# Create scatter plot of BETA vs BETA.AFR, colored by P.AFR threshold
afrPlot <- ggplot(allSentsBonf, aes(x = BETA.EUR, y = BETA.AFR, color = P.AFR < thresh)) +
  geom_point() + 
  geom_errorbar(aes(ymin = BETA.AFR - SE.AFR , ymax = BETA.AFR + SE.AFR)) + 
  geom_errorbarh(aes(xmin = BETA.EUR - SE.EUR, xmax = BETA.EUR + SE.EUR))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
  scale_color_manual(values = c("TRUE" = "#CC0066", "FALSE" = "#6699CC")) +
  labs(x = "BETA (EUR)", y = "BETA (AFR)", color =  substitute(paste("P"[AFR], "<", threshLabel), list(threshLabel = threshLabel)))

# Create scatter plot of BETA vs BETA.CSA, colored by P.CSA threshold
csaPlot <- ggplot(allSentsBonf, aes(x = BETA.EUR, y = BETA.CSA, color = P.CSA < thresh)) +
  geom_point() +
  geom_errorbar(aes(ymin = BETA.CSA - SE.CSA , ymax = BETA.CSA + SE.CSA)) + 
  geom_errorbarh(aes(xmin = BETA.EUR - SE.EUR, xmax = BETA.EUR + SE.EUR))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
  scale_color_manual(values = c("TRUE" = "#CC0066", "FALSE" = "#6699CC")) +
  labs(x = "BETA (EUR)", y = "BETA (CSA)", color =  substitute(paste("P"[CSA], "<", threshLabel), list(threshLabel = threshLabel)))

  png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/plots/ancestriesComparisonBonfSigFPCSentinels.png", width = 800, height = 400)
  afrPlot / csaPlot %>%
  print
  dev.off()


fwrite(allSentsBonf, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/ancestryComparisonsFPCs.csv")


allSentsBonf[, EurAfrConcordant] %>% table
allSentsBonf[, EurCsaConcordant] %>% table

allSentsBonf[, .(FPC, EurAfrConcordant)] %>% table
allSentsBonf[, .(FPC, EurCsaConcordant)] %>% table

cor.test(allSentsBonf[, BETA.EUR],  allSentsBonf[, BETA.AFR])
cor.test(allSentsBonf[, BETA.EUR],  allSentsBonf[, BETA.CSA])
