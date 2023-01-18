#!/usr/bin/env Rscript

library(data.table)
library(magrittr)
library(tidyverse)
library(corrplot)

chr <- 22

pixels <- fread("/vast/scratch/users/jackson.v/retThickness/GWAS/pixels.txt") %>%
setnames(., c("pixel", "y", "x"))

## full sig results
dir <- paste0("/vast/scratch/users/jackson.v/retThickness/GWAS/GWsigResults/chr",chr)

results <- lapply(c(1:119), function(slice) {

  print(paste(slice))

  file <- paste0(dir,"/chr",chr,"Slice",slice,"_5e-5Sig.txt")

  sliceResults <- fread(file)
  return(sliceResults)

}) %>%
rbindlist


## grid - n snps for each pixel
dt <- results[P<5e-8] %>%
  .[ , .N, by = pixel] %>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert = T)] %>%
  .[!is.na(y)]

plot <- ggplot(dt) +
    geom_tile(aes(x = x, y = y, fill = N)) +
    scale_fill_gradient2() +
    scale_y_reverse() +
    theme_bw() +
    theme(legend.position = "bottom")

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/plots/fullGrid_chr",chr,"_pixelGWSigSNPs.png"), width = 1200, height = 600)
print(plot)
dev.off()


# ## N spixels for each SNP
# dt <- results[ , .N, by = POS] %>%
#   .[,POS := as.numeric(POS)]

# png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/plots/fullGrid_chr",chr,"_ChromosomeSigSNPs.png"), width = 1200, height = 600)
# ggplot(dt, aes(x=POS, y = N)) +
# geom_line()
# dev.off()

## full clumped resultsDir  clumped <- lapply(pixels[, pixel], function(pix) {
dir <- "/vast/scratch/users/jackson.v/retThickness/GWAS/clumpedResults/chr22/"
clumped <- lapply(pixels[, pixel], function(pix) {
    
    print(paste0(pix))

    slice <- pixels[pixel == pix, y]

    file <- paste0(dir,"/",slice,"/chr",chr,"Pixel.",pix,".clumped")

    if(file.exists(file)) {
      pixResults <- fread(file, fill=T) %>%
      .[, pixel := pix] %>%
      .[!is.na(CHR)]

    if(nrow(pixResults) > 0) {
      return(pixResults)
    }

    }

  }) %>%
  rbindlist





results.filt <- results[P<5e-8]
sentinels <- NULL
allocatedCheck <- NULL


while(nrow(results.filt)>0) {
		# select top hit
		hit <- results.filt[P==min(results.filt[,P])]

    print(paste(hit[,ID], " with P =", hit[,P]))

    snp <- hit[1,ID]
    pix <- hit[1,pixel]
    clumpedResult <- clumped[SNP==snp & pixel==pix]

    all.pixels <- results.filt[ID==snp, pixel] %>% unique

    if(nrow(clumpedResult)==0) {
      print(paste("SNP", snp, "Not sentinel for chunk"))
      cat("\n")
    
      clumpedResult2 <- clumped[SP2 %like% snp & pixel==pix]
      clumpSentinel <- clumpedResult2[,SNP]

     if(nrow(sentinels[SNPsInLocus %like% clumpSentinel & pixels %like% pix])>0) {

      print(paste("SNP", snp, "already assigned to locus"))
      cat("\n")
     } else {

      check <- data.table(snp=snp,
        pixels= paste(all.pixels, collapse = ","),
        nPix = length(all.pixels))

      allocatedCheck <- rbind(allocatedCheck, check)

     }

    results.filt <- results.filt[!(ID %in% snp & pixel %in% all.pixels)]

    } else{

    tmp <- clumpedResult[,SP2] %>%
    strsplit(., ",") %>%
    unlist
    clumpSNPs <- sapply(strsplit(tmp,"\\("), `[`, 1)

    clumpPixels <- results.filt[ID %in% c(snp,clumpSNPs), pixel] %>% unique
    n <- length(clumpSNPs)
    noSNPs <- clumpSNPs[1]=="NONE"
    nSNPs <- ifelse(noSNPs, 0, n)
    nPix <- length(clumpPixels)

		hit <- hit[, nSNPsLocus := nSNPs] %>%
      .[, nPixelsLocus := nPix] %>%
      .[, SNPsInLocus := paste(clumpSNPs, collapse = ",")] %>%
     .[, pixels := paste(clumpPixels, collapse = ",")]

		# remove SNPs in clump
		results.filt <- results.filt[!(ID %in% c(snp,clumpSNPs) & pixel %in% c(pix,clumpPixels))]

		# Add hit onto list of sentinel SNPs
		sentinels <- rbind(sentinels, hit[1,])

    print(paste(nSNPs, "additional SNPs in locus"))
    print(paste(nPix, "additional Pixels significant at locus"))
    print(paste(nrow(results.filt), "SNPs remaining..."))
    cat("\n")
}
}


fwrite(sentinels, file ="/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/chr22sentinels.csv", sep = ",")






mat <- results[ID %in% sentinels[,ID]] %>%
dcast(., pixel ~ ID, value.var = "BETA", fun.aggregate=mean) %>%
  as.matrix(., rownames = "pixel")
mat[is.na(mat)] <- 0
class(mat) <- "numeric"

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/plots/chr",chr,"_senitnelsBetasCorrPlot.png"), width = 1200, height = 1200)
mat %>% cor %>% corrplot
dev.off()

