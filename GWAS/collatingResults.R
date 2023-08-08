#!/usr/bin/env Rscript

library(data.table)
library(magrittr)
library(tidyverse)
library(corrplot)
library(patchwork)
library(optparse)

option_list <-  list(
  make_option(c("-c", "--chr"), type="integer", default=NULL,
              help="chromosome number", metavar="integer"));

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

chr <-  ifelse(opt$chr==23, "X", opt$chr)


pixels <- fread("/vast/scratch/users/jackson.v/retThickness/GWAS/pixels.txt") %>%
setnames(., c("pixel", "y", "x"))

## full sig results
dir <- paste0("/vast/scratch/users/jackson.v/retThickness/GWAS/GWsigResults/chr",chr)

results <- lapply(c(1:119), function(slice) {

  print(paste(slice))

  file <- paste0(dir,"/chr",chr,"Slice",slice,"_5e-5Sig.txt")

  sliceResults <- fread(file) %>%
    setnames(., "#POS", "POS")
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


# ## N pixels for each SNP
dt <- results[ , .N, by = POS] %>%
  .[,POS := as.numeric(POS)]

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/plots/fullGrid_chr",chr,"_ChromosomeSigSNPs.png"), width = 1200, height = 600)
ggplot(dt, aes(x=POS, y = N)) +
geom_line()
dev.off()

## full clumped resultsDir  clumped <- lapply(pixels[, pixel], function(pix) {

  print("read in clumped results")
dir <- paste0("/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/clumpedResults/chr",chr)
clumped <- lapply(pixels[, pixel], function(pix) {
    
    # print(paste0(pix))

    slice <- pixels[pixel == pix, y]

    file <- paste0(dir,"/",slice,"/chr",chr,"Pixel.",pix,"_thresh0.001_withOverlap.clumped")

    if(file.exists(file)) {
      print(paste0(pix))
      pixResults <- fread(file, fill=T) %>%
      .[, pixel := pix] %>%
      .[!is.na(CHR)]

    if(nrow(pixResults) > 0) {
      return(pixResults)
    }

    }

  }) %>%
  rbindlist





results.filt.full <- results[P<5e-5]
results.filt <- results.filt.full

sentinels <- NULL
allocatedCheck <- NULL


while(nrow(results.filt[P<5e-8])>0) {
		# select top hit
		hit <- results.filt[P==min(results.filt[,P])]

    print(paste(hit[,ID], " with P =", hit[,P]))

    snp <- hit[1,ID]
    pix <- hit[1,pixel]
    clumpedResult <- clumped[SNP==snp & pixel==pix]

    all.pixels <- results.filt.full[ID==snp, pixel] %>% unique

    if(nrow(clumpedResult)==0) {
      print(paste("SNP", snp, "Not sentinel for chunk"))
      cat("\n")

      # find clump
      pixClumps  <- clumped[pixel==pix]

      whichClump <- lapply(pixClumps[,SNP], function(sent) {

          tmp <- pixClumps[SNP==sent,SP2] %>%
            strsplit(., ",") %>%
            unlist
          clumpSNPs <- sapply(strsplit(tmp,"\\("), `[`, 1)

        inClump <-  ifelse(snp %in% clumpSNPs, T, F)
        return(inClump)

      }) %>%
      unlist 

      clumpedResult2 <- pixClumps[whichClump]
      clumpSentinel <- clumpedResult2[,SNP]
       tmp <- pixClumps[whichClump,SP2] %>%
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
      print(paste("removing other clump SNPs;", nrow(results.filt[!(ID %in% clumpSNPs & pixel %in% all.pixels) & P<5E-8]), "SNPs remaining..."))
      cat("\n")

      check <- data.table(snp=snp,
        clumpTopSNPPOS=clumpedResult2[,BP],
        clumptopSNP=clumpSentinel,
        topPixel=pix,
        nPix = length(all.pixels),
        assignedSentinelPOS = sentinels[whichSentinel, POS],
        assignedSentinelSNP = sentinels[whichSentinel, ID],
        assignedSentinelPixel = sentinels[whichSentinel, pixel],
        allocated = "Y",
        pixels= paste(all.pixels, collapse = ","),
        otherClumpSNPs=clumpedResult2[,SP2])

      allocatedCheck <- rbind(allocatedCheck, check)


     } else {

      check <- data.table(snp=snp,
        clumpTopSNPPOS=clumpedResult2[,BP],
        clumptopSNP=clumpSentinel,
        topPixel=pix,
        nPix = length(all.pixels),
        assignedSentinelPOS = NA,
        assignedSentinelSNP = NA,
        assignedSentinelPixel = NA,
        allocated = "N",
        pixels= paste(all.pixels, collapse = ","),
        otherClumpSNPs=clumpedResult2[,SP2])

      allocatedCheck <- rbind(allocatedCheck, check)

     }

    results.filt <- results.filt[!(ID %in% clumpSNPs & pixel %in% all.pixels)]

    } else {

    tmp <- clumpedResult[,SP2] %>%
    strsplit(., ",") %>%
    unlist
    clumpSNPs <- sapply(strsplit(tmp,"\\("), `[`, 1)

    clumpPixels <- results.filt.full[ID %in% c(snp,clumpSNPs), pixel] %>% unique
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
    print(paste(nrow(results.filt[P<5E-8]), "SNPs remaining..."))
    cat("\n")
}
}


fwrite(sentinels, file = paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/chr",chr,"sentinels_clumpThresh0.001_withOverlap.txt"), sep = "\t")
fwrite(sentinels[,.(ID)], file = paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/chr",chr,"sentinelsIDonly_clumpThresh0.001_withOverlap.txt"), col.names=F)

if(!is.null(allocatedCheck)) {
  fwrite(allocatedCheck, file = paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/chr",chr,"_checkAllocated_clumpThresh0.001_withOverlap.txt"), sep = "\t")
}

order <- sentinels[,.(ID, POS)] %>%
unique %>%
setorder(., POS)  %>%
.[, ID]

mat <- results[ID %in% sentinels[,ID]] %>%
dcast(., pixel~ ID , value.var = "BETA", fun.aggregate=mean) %>%
  as.matrix(., rownames = "pixel")
mat[is.na(mat)] <- 0
class(mat) <- "numeric"

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/plots/chr",chr,"/chr",chr,"_sentinelsBetasCorrPlot_clumpThresh0.001_withOverlap.png"), width = 1200, height = 1200)
mat[,order] %>% cor %>% corrplot(., order = "original")
dev.off()


sentinelIDs <- sentinels[,ID] %>% unique

lapply(sentinelIDs, function(snp) {

snpOut <- ifelse(snp %like% ":", str_replace(snp, ":", "_"), snp)

snpResult <- results[ID == snp] %>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert = T)] %>%
  .[!is.na(y)] %>%
  .[, log10P := log(P, 10)]


statsPlots <- lapply(c("BETA", "T_STAT", "log10P"), function(stat) {
    plot <- ggplot(snpResult) +
    geom_tile(aes_string(x = "x", y = "y", fill = stat)) +
    # geom_path(aes(x=col,y=row,color=area),data = areas, size = 0.5) +
    scale_fill_gradient2() +
    scale_y_reverse() +
    theme_bw() +
    theme(legend.position = "bottom")+
    ggtitle(paste(snp,stat, snp <-sep = " - "))

  return(plot)
})

  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/plots/chr",chr,"/assocsPixelwise_",snpOut,".png"), width = 1800, height = 600)
  reduce(statsPlots , `+`) %>%
  print
  dev.off()

})

