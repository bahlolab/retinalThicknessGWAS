#!/usr/bin/env Rscript

library(data.table)
library(magrittr)
library(tidyverse)
library(here)
library(patchwork)
library(optparse)

option_list <-  list(
  make_option(c("-c", "--chr"), type="integer", default=NULL,
              help="chromosome number", metavar="integer")
);

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

chr <-  opt$chr


dataDir <- paste0("/vast/scratch/users/jackson.v/retThickness/GWAS/sentinelResults/chr",chr)

results <- lapply(c(1:119), function(slice) {
    paste0(dataDir,"/chr",chr,"Slice",slice,"_sentinels.txt") %>%
    fread
  }) %>%
  rbindlist %>%
  # setnames(., c("CHROM", "POS", "ID", "REF", "ALT", "A1", "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P", "ERRCODE", "pixel")) %>%
  .[, log10P := (-1)*log10(P)] %>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert=TRUE)]


sentinels <- results[,ID] %>% unique

lapply(sentinels, function(snp) {

snpResult <- results[ID == snp]

statsPlots <- lapply(c("BETA", "T_STAT", "log10P"), function(stat) {
    plot <- ggplot(snpResult) +
    geom_tile(aes_string(x = "x", y = "y", fill = stat)) +
    # geom_path(aes(x=col,y=row,color=area),data = areas, size = 0.5) +
    scale_fill_gradient2() +
    scale_y_reverse() +
    theme_bw() +
    theme(legend.position = "bottom")+
    ggtitle(paste(snp,stat, sep = " - "))

  return(plot)
})
  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/plots/chr",chr,"/chr",chr,"_",snp,"AssocsPixelwise.png"), width = 1800, height = 600)
  reduce(statsPlots , `+`) %>%
  print
  dev.off()

})





























## misc old code

# #snpResult <- lapply(c(1:nrow(pixels)), function(idx) {
# snpResult <- lapply(c(1:10), function(idx) {
#
#   slice <- pixels[idx, y]
#   pixel <- pixels[idx, pixel]
#
#   print(paste(pixel))
#
#  result <- fread(paste0("test/",slice,"/pixel.",pixel,".glm.linear")) %>%
#   .[TEST == "ADD"] %>%
#   .[ID %in% snp] %>%
#   .[, pixel := pixel]
#
#   return(result)
#
# }) %>%
# rbindlist
#
#
#
# ## single results file
# snpResult <- fread("/vast/scratch/users/jackson.v/retThickness/GWAS/test/rs17421627Result.txt") %>%
#   .[, log10P := (-1)*log10(P)] %>%
#   .[, c("ERRCODE", "pixel") := tstrsplit(get("ERRCODE pixel"), " ")] %>%
#   .[, c("y", "x") := tstrsplit(pixel, "_", type.convert=TRUE)]
#
#
# for(stat in c("BETA", "T_STAT", "log10P")) {
#
#   plot <- ggplot(snpResult , aes_string(x = "x", y = "y", fill = stat)) +
#     geom_tile() +
#     scale_fill_gradient2()
#
#   png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/testPlots/rs1742162AssocsPixelwise_",stat,".png"))
#   print(plot)
#   dev.off()
#
# }
#
#
#
