#!/usr/bin/env Rscript

library(data.table)
library(magrittr)
library(tidyverse)
library(optparse)

option_list <-  list(
  make_option(c("-c", "--chr"), type="integer", default=NULL,
              help="chromosome number", metavar="integer"),
  make_option(c("-s", "--snp"), type="character", default=NULL,
              help="SNP", metavar="character")
  
);

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

chr <-  ifelse(opt$chr==23, "X", opt$chr)
SNP <- opt$snp

dir <- paste0("/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/results/chr",chr)
sentinels <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/chr",chr,"sentinels_clumpThresh0.001_withOverlap.txt") %>%
 fread 


    pix <- sentinels[ID==SNP, pixel]
    slice <- str_split(pix, "_") %>%
        unlist %>% 
        .[1]

    pos <- sentinels[ID==SNP, POS]
    posMin <- pos - 1000000
    posMax <- pos + 1000000

    results <- paste0(dir,"/",slice,"/chr",chr,"Pixel.",pix,".glm.linear.gz") %>% 
    fread %>%
    setnames(., "#POS", "POS") %>%
    .[POS %between% c(posMin, posMax)] %>%
    .[, .(ID, P)]

    # locuszoom.commands <- data.table(snp=SNP,
    #                              chr=NA,
    #                              start=NA,
    #                              end=NA,
    #                              flank="500kb",
    #                              run="yes",
    #                              arguments="showAnnot=T showRecomb=T")


    fwrite(results, file = paste0("/vast/scratch/users/jackson.v/retThickness/GWAS/regionPlots/chr",chr,"/",SNP,"_METAL.txt"), sep = "\t")
    # fwrite(locuszoom.commands, file = "", sep = "\t")


