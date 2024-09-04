#!/usr/bin/env Rscript


library(data.table)
library(magrittr)
library(ggplot2)
library(stringr)
library(optparse)

option_list <-  list(
  make_option(c("-c", "--chr"), type="integer", default=NULL,
              help="chromosome number", metavar="integer"),
  make_option(c("-d", "--directory"), type="character", default=NULL,
              help="clumped results directory", metavar="character")
);

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

chr <-  opt$chr
dir <- opt$directory

#chr <- 22
#dir <- "/vast/scratch/users/jackson.v/retThickness/GWAS/clumpedResults"

## clumped results to be used for initial filtering 
## clumps of SNPs identified as those with P < 5e-5).

## No filtering at this stage.

## older versions included the following:
## identify signals with few supporting SNPs (<5 n clump).
## identify pixels with an excess of (clumped) signals
## identify pixels with excess of signals with few supporting SNPs.

pixels <- fread("/vast/scratch/users/jackson.v/retThickness/GWAS/pixels.txt") %>%
  setnames(., c("pixel", "y", "x"))


clumped <- lapply(pixels[, pixel], function(pix) {
  
  slice <- pixels[pixel == pix, y]
  file <- paste0(dir,"/chr",chr,"/",slice,"/chr",chr,"Pixel.",pix,"_thresh0.001_withOverlap.clumped")
  
  if(file.exists(file)) {
    pixResults <- fread(file, fill=T) 
    
    
    topSNP <- pixResults[,SNP]
    clumpSNPs <- pixResults[,SP2] 
    
    nSig5e5 <-  sapply(clumpSNPs, function(S) { 
      strsplit(S, ",") %>% 
        unlist %>% 
        length 
    })
    
    nSig1e4 <- pixResults[,S0001] 
    
    pixResults <- data.table(topSNP, clumpSNPs, nSig5e5, nSig1e4) %>%
      .[, pixel := pix] %>%
      .[!is.na(nSig1e4)]
    
    
    if(length(pixResults) > 0) {
      return(pixResults)
    }
    
  }
  
}) %>%
  rbindlist  %>%
  .[, chr := chr]


## summary of nuber of clumps per pixel
pixClumped <- clumped %>%
  .[ , .N, by = pixel]

p1 <- ggplot(pixClumped, aes(x = N)) +
  geom_histogram(bins = 50) +
  scale_y_continuous(trans='log10')

# png(paste0("stornext/Bioinf/data/lab_bahlo/projects/misc/retinalThickness/GWAS/output/pixelWiseResultsPlots/chr",chr,"_clumpsPerPixel.png"))
# print(p1)
# dev.off()

## number of "singleton SNPs" per pixel
singletons <- clumped %>%
  .[, sum(nSig1e4==0, na.rm = T) , by = pixel]

p2 <- ggplot(singletons, aes(x = V1)) +
  geom_histogram()+
  scale_y_continuous(trans='log10')

# png(paste0("stornext/Bioinf/data/lab_bahlo/projects/misc/retinalThickness/GWAS/output/pixelWiseResultsPlots/chr",chr,"_singeltons1e-4PerPixel.png"))
# print(p2)
# dev.off()

# number of SNPs per clumps
p3 <- ggplot(clumped[nSig1e4 > 0], aes(x = nSig1e4)) +
  geom_histogram(bins = 50) +
  scale_x_continuous(trans='log10')

# png(paste0("stornext/Bioinf/data/lab_bahlo/projects/misc/retinalThickness/GWAS/output/pixelWiseResultsPlots/chr",chr,"_snps1e-4PerClump.png"))
# print(p3)
# dev.off()

p4 <- ggplot(clumped[clumpSNPs != "NONE"], aes(x = nSig5e5)) +
  geom_histogram(bins = 50) +
  scale_x_continuous(trans='log10')

# png(paste0("stornext/Bioinf/data/lab_bahlo/projects/misc/retinalThickness/GWAS/output/pixelWiseResultsPlots/chr",chr,"_snps5e-5PerClump.png"))
# print(p4)
# dev.off()



tmp <- clumped[,clumpSNPs] %>%
  strsplit(., ",") %>%
  unlist %>%
  unique

clumpSNPs <- sapply(strsplit(tmp,"\\("), `[`, 1)

uniqueSNPs <- c(clumped[,topSNP], clumpSNPs) %>% unique
fwrite(data.table(ID=uniqueSNPs), file = paste0("/vast/scratch/users/jackson.v/retThickness/GWAS/GWsigResults/chr",chr,"/chr",chr,"all5e-5SigSNPs.txt"), row.names=F, col.names=F)


