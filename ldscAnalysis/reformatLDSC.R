#!/usr/bin/env Rscript

library(optparse)
library(data.table)
library(magrittr)
library(here)

# Parse command line arguments
option_list <-  list(
 make_option(c("-s", "--slice"), type="integer", default=NULL,
              help="slice number", metavar="integer"));

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

slice <- opt$slice

pixels <- fread(here("pixels.txt")) %>%
  setnames(., c("pixel", "y", "x")) %>%
  .[y==slice]


refs <- lapply(c(1:22), function(chr) {
  here(paste0("chr",chr,"SNPinfo.txt")) %>%
  fread %>%
  setnames(., c("SNP", "CHR", "BP", "A1", "A2", "freq", "GENE"))
})

lapply(pixels[,pixel], function(pix) {

outName <- paste0(pix, "_ldscFormat.txt") %>%
here::here("ldscInFiles", slice,  .)

pixResult <- lapply(c(1:22), function(chr) {

ref <- refs[[chr]]

data <- here::here("results", paste0("chr",chr), slice, paste0("chr",chr,"Pixel.",pix,".glm.linear.gz")) %>%
  fread(.) %>%
  setnames(., c("BP", "SNP", "A1", "b", "se", "Tstat", "p")) %>%
  .[ref, on = c("SNP", "A1")] %>%
  .[, N := 43148] %>%
  .[, .(SNP, A1, A2, N, freq, p, b)] %>%
  setnames(., "freq", "FRQ") %>%
  setnames(., "b", "BETA") %>%
  setnames(., "p", "P")

  return(data)

}) %>%
rbindlist

print(paste0("Writing to ", outName))
fwrite(pixResult, outName, sep = " ")

rm(pixResult)

})




