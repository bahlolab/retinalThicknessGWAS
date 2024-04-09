library(data.table)
library(magrittr)


pixels <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/pixels.txt") %>%
setnames(., c("pixel", "y", "x"))


## full sig results
dir <- "/vast/scratch/users/jackson.v/retThickness/AFRGWAS/results"

results <- lapply(c(1:22, "X"), function(chr) {

print(paste(chr))
chrResult <- lapply(c(1:119), function(slice) {

#   print(paste(slice))

pix <- pixels[y==slice, pixel]


  sliceResults <- lapply(pix, function(i) {
  
    pixResult <- paste0(dir,"/chr",chr,"/",slice,"/chr",chr,"PixelAFR.",i,".glm.linear") %>%
        fread(.) %>%
        setnames(., "#CHROM", "CHROM") %>%
        .[, .(CHROM, POS, ID, A1, A1_FREQ, BETA, SE, P)] %>%
        .[, pixel := i]

    return(pixResult)

  }) %>%
  rbindlist

  return(sliceResults)

}) %>%
rbindlist

fwrite(chrResult, file = paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinelsAFR/chr",chr,"EURsentinels_AFRassocsAllPixels.txt"))

})


dir <- "/vast/scratch/users/jackson.v/retThickness/AFRGWAS/results"

results <- lapply(c(1:22, "X"), function(chr) {

print(paste(chr))
chrResult <- lapply(c(1:119), function(slice) {

#   print(paste(slice))

pix <- pixels[y==slice, pixel]


  sliceResults <- lapply(pix, function(i) {
  
    pixResult <- paste0(dir,"/chr",chr,"/",slice,"/chr",chr,"PixelCSA.",i,".glm.linear") %>%
        fread(.) %>%
        setnames(., "#CHROM", "CHROM") %>%
        .[, .(CHROM, POS, ID, A1, A1_FREQ, BETA, SE, P)] %>%
        .[, pixel := i]

    return(pixResult)

  }) %>%
  rbindlist

  return(sliceResults)

}) %>%
rbindlist

fwrite(chrResult, file = paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinelsCSA/chr",chr,"EURsentinels_CSAassocsAllPixels.txt"))

})


