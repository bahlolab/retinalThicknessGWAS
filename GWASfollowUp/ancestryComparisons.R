library(data.table)
library(magrittr)
library(tidyverse)
library(patchwork)

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
.[A1_FREQ < 0.001, BETA := NA]


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
.[A1_FREQ < 0.001, BETA := NA]


allSents <- sents %>%
    merge(., sentsAFR, by = c("CHROM", "POS", "ID", "A1", "pixel"), suffixes = c("", ".AFR")) %>%
    merge(., sentsCSA, on = c("POS", "ID", "REF", "ALT", "A1", "pixel"), suffixes = c("", ".CSA"))


allSentsBonf <- allSents[BonferroniSig == "Y"]
thresh <- 0.05/nrow(allSentsBonf)

# Create scatter plot of BETA vs BETA.AFR, colored by P.AFR threshold
afrPlot <- ggplot(allSentsBonf, aes(x = BETA, y = BETA.AFR, color = P.AFR < thresh)) +
  geom_point() +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  labs(x = "BETA", y = "BETA.AFR", color = "P.AFR < thresh")

# Create scatter plot of BETA vs BETA.CSA, colored by P.CSA threshold
csaPlot <- ggplot(allSentsBonf, aes(x = BETA, y = BETA.CSA, color = P.CSA < thresh)) +
  geom_point() +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  labs(x = "BETA", y = "BETA.CSA", color = "P.CSA < thresh")

  png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/plots/ancestriesComparisonBonfSigSentinels.png", width = 1200, height = 600)
  afrPlot + csaPlot %>%
  print
  dev.off()


thresh <- 0.05/nrow(allSents)

# Create scatter plot of BETA vs BETA.AFR, colored by P.AFR threshold
afrPlot <- ggplot(allSents, aes(x = BETA, y = BETA.AFR, color = P.AFR < thresh)) +
  geom_point() +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  labs(x = "BETA", y = "BETA.AFR", color = "P.AFR < thresh")

# Create scatter plot of BETA vs BETA.CSA, colored by P.CSA threshold
csaPlot <- ggplot(allSents, aes(x = BETA, y = BETA.CSA, color = P.CSA < thresh)) +
  geom_point() +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  labs(x = "BETA", y = "BETA.CSA", color = "P.CSA < thresh")

  png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/plots/ancestriesComparisonGWSigSentinels.png", width = 1200, height = 600)
  afrPlot + csaPlot %>%
  print
  dev.off()