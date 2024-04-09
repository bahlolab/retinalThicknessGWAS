#!/usr/bin/env Rscript

library(optparse)
library(data.table)
library(magrittr)
library(here)

# Parse command line arguments
option_list <-  list(
  make_option(c("-c", "--chr"), type="integer", default=NULL,
              help="chromosome number", metavar="integer"),
 make_option(c("-s", "--slice"), type="integer", default=NULL,
              help="slice number", metavar="integer"),
 make_option(c("-p", "--pixel"), type="character", default=NULL,
              help="pixel", metavar="character"));

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

chr <-  opt$chr
slice <- opt$slice
pixel <- opt$pixel

resultsChr <- ifelse(chr == 23, "X", chr)
file <-  here::here("results", paste0("chr",resultsChr), slice, paste0("chr",resultsChr,"Pixel.",pixel,".glm.linear.gz")) 

outName <- paste0("chr",chr,"Pixel.",pixel, "_cojoFormat.txt") %>%
here::here("cojoGWsigInFiles", .)

ref <- here(paste0("chr",resultsChr,"SNPinfo.txt")) %>%
fread %>%
setnames(., c("SNP", "CHR", "BP", "A1", "A2", "freq", "GENE")) %>%
.[, CHR := chr]


  print(paste0("Processing ", file))
  data <- fread(file) %>%
  setnames(., c("BP", "SNP", "A1", "b", "se", "Tstat", "p")) %>%
  .[ref, on = c("SNP", "A1")] %>%
  .[, N := 43148] %>%
  .[, .(SNP, A1, A2, freq, b, se, p, N)]

  print(paste0("Writing to ", outName))
  fwrite(data, outName, sep = " ")
