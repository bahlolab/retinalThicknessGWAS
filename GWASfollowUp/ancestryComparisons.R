library(data.table)
library(magrittr)
library(tidyverse)
library(patchwork)
library(scales)

pixels <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/pixels.txt") %>%
setnames(., c("pixel", "y", "x"))

sents <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/annotatedSentinelsPixelwise.txt") %>%
    .[, CHROM := ifelse(CHR==23, "X", CHR)] %>%
    .[, .(CHROM, POS, ID, EffAllele, EffAlleleFreq, BETA, P, BonferroniSig, topPixel)] %>%
    setnames(., c("CHROM", "POS", "ID", "A1", "A1_FREQ", "BETA", "P","BonferroniSig", "pixel" ))

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
    .[, c("CHROM", "POS", "ID", "A1", "A1_FREQ", "BETA", "P", "pixel" )]

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
    .[, c("CHROM", "POS", "ID", "A1", "A1_FREQ", "BETA", "P", "pixel" )]

    return(chrResult)}) %>%

rbindlist %>%
.[A1_FREQ < 0.001, c("BETA", "P") := NA]

allSents <- sents %>%
    merge(., sentsAFR, by = c("CHROM", "POS", "ID", "pixel"), suffixes = c("", ".AFR"), all = T) %>%
    merge(., sentsCSA, on = c("POS", "ID", "REF", "ALT", "pixel"), suffixes = c("", ".CSA"), all =T) %>%
    setnames(., old = c("A1_FREQ", "BETA", "P"), new = c("A1_FREQ.EUR", "BETA.EUR", "P.EUR")) %>%
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
   .[, .(CHROM, POS, ID, pixel, A1, A1_FREQ.EUR, BETA.EUR, P.EUR, A1_FREQ.AFR, BETA.AFR, P.AFR, A1_FREQ.CSA, BETA.CSA, P.CSA, EurAfrConcordant, EurCsaConcordant)]          

thresh <- 0.05/nrow(allSentsBonf)
threshLabel <- scientific(thresh, digits = 3)

# Create scatter plot of BETA vs BETA.AFR, colored by P.AFR threshold
afrPlot <- ggplot(allSentsBonf, aes(x = BETA.EUR, y = BETA.AFR, color = P.AFR < thresh)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
  scale_color_manual(values = c("TRUE" = "#CC0066", "FALSE" = "#6699CC")) +
abs(x = "BETA (EUR)", y = "BETA (AFR)", color =  substitute(paste("P"[AFR], "<", threshLabel), list(threshLabel = threshLabel)))

# Create scatter plot of BETA vs BETA.CSA, colored by P.CSA threshold
csaPlot <- ggplot(allSentsBonf, aes(x = BETA.EUR, y = BETA.CSA, color = P.CSA < thresh)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
  scale_color_manual(values = c("TRUE" = "#CC0066", "FALSE" = "#6699CC")) +
labs(x = "BETA (EUR)", y = "BETA (CSA)", color =  substitute(paste("P"[CSA], "<", threshLabel), list(threshLabel = threshLabel)))

  png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/plots/ancestriesComparisonBonfSigSentinels.png", width = 1200, height = 600)
  afrPlot + csaPlot %>%
  print
  dev.off()


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
