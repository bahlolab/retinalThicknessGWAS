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
              help="slice number", metavar="integer"));

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

chr <-  opt$chr
slice <- opt$slice

files <- list.files(path = here::here("results", paste0("chr",chr), slice), 
  pattern = "*.glm.linear.gz",
  full.names = TRUE) 

outNames <- gsub(".*/(chr\\d+Pixel\\.\\d+_\\d+)\\.glm\\.linear\\.gz", "\\1", files) %>%
paste0(., "_cojoFormat.txt") %>%
here::here("cojoInFiles", paste0("chr",chr), slice,  .)

ref <- here(paste0("chr",chr,"SNPinfo.txt")) %>%
fread %>%
setnames(., c("SNP", "CHR", "BP", "A1", "A2", "freq", "GENE"))

for (i in 1:length(files)) {
  file <- files[i]
  outName <- outNames[i]
  print(paste0("Processing ", file))
  data <- fread(file) %>%
  setnames(., c("BP", "SNP", "A1", "b", "se", "Tstat", "p")) %>%
  .[ref, on = c("SNP", "A1")] %>%
  .[, N := 43148] %>%
  .[, .(SNP, A1, A2, freq, b, se, p, N)]

  print(paste0("Writing to ", outName))
  fwrite(data, outName, sep = " ")
}
