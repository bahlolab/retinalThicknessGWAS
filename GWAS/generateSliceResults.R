library(data.table)
library(magrittr)
library(ggplot2)
library(here)
library(optparse)

option_list <-  list(
  make_option(c("-c", "--chr"), type="integer", default=NULL,
              help="chromosome, integer 1-22", metavar="integer"),
  make_option(c("-s", "--slice"), type="integer", default=NULL,
              help="slice number, integer 1-119", metavar="integer")
);

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

chr <- opt$chr
slice <- opt$slice


pixels <- fread("/vast/scratch/users/jackson.v/retThickness/GWAS/pixels.txt") %>%
  setnames(., c("pixel", "y", "x")) %>%
  .[y == slice]

pix <- pixels[,pixel]

pixResults <- lapply(pix, function(pixel) {


  print(paste(pixel))

  file <- paste0("/vast/scratch/users/jackson.v/retThickness/GWAS/results/chr",chr,"/",slice,"/chr",chr,"Pixel.",pixel,".glm.linear.gz")
 # if(file.exists(file)) {

    result <- fread(file, select = c(1:4, 7)) %>%
      setnames(., c("POS", "ID", "A1", paste0(pixel,"_BETA"), paste0(pixel,"_P"))) %>%
      setkey(., POS, ID, A1)

    return(result)

#  }
}) 

myMerge <- function(x, y) {
  x[y]
}

sliceResult <- Reduce(myMerge, pixResults)

fwrite(sliceResult, file = paste0("/vast/scratch/users/jackson.v/retThickness/GWAS/results/chr",chr,"/",slice,"_result.txt.gz"), compress =  "gzip")

