library(data.table)
library(magrittr)
library(tidyverse)



chroms <- c(1:22)
fpcs <- c(1:50)

####################
## rough code

## read in results and write out FUMA input file
lapply(fpcs, function(i) {
  
result <- lapply(chroms, function(chr) {
  
res <- paste0("/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/fpcGWAS/results/chr",chr,"/chr",chr,"EUR.fpc",i,".glm.linear") %>%
  fread %>%
  setnames(., "#CHROM", "CHR") %>%
#  .[A1_CT >= 500]
return(res)

}) %>%
rbindlist



FUMAResults <- result %>%
  .[, .(CHR, POS, REF, ALT, OBS_CT, BETA, P)]
fwrite(FUMAResults, file =paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWAS/output/FUMAresults/fpc",i,"_FUMA.txt"), sep="\t")

})




i <- 1

result <- lapply(chroms, function(chr) {
  
  res <- paste0("/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/fpcGWAS/results/chr",chr,"/chr",chr,"EUR.fpc",i,".glm.linear") %>%
    fread %>%
    setnames(., "#CHROM", "CHR") %>%
    setnames(., "POS", "BP") %>%
    setnames(., "ID", "SNP") %>%
    setnames(., "REF", "A2") %>%
    setnames("A1_FREQ", "FREQ") %>%
    setnames(., "OBS_CT", "N") %>%
    return(res)
  
}) %>%
  rbindlist



outResults <- result %>%
  .[, .(CHR, BP, SNP, A1, A2, FREQ, BETA, SE, P, N)]
fwrite(outResults, file =paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWAS/output/FUMAresults/fpc",i,"_CTG-VLtest.txt"), sep="\t")
