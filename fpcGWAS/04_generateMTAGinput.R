library(data.table)
library(magrittr)
library(tidyverse)



chroms <- c(1:22)
fpcs <- c(1:50)
snpAlleles <- c("A","C", "G", "T")

####################
## rough code

## read in results and write out FUMA input file
lapply(fpcs, function(i) {
  
print(paste(i))

result <- lapply(chroms, function(chr) {

res <- paste0("/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/fpcGWAS/results/chr",chr,"/chr",chr,".fpc",i,".glm.linear") %>%
  fread %>%
  setnames(., "#CHROM", "CHR") %>%
  .[A1_CT >= 500] %>%
  .[,A2 := ifelse(A1==REF, ALT, REF)] %>%
  .[A1 %in% snpAlleles & A2 %in% snpAlleles] 

  return(res)

}) %>%
rbindlist

MTAGResults <- result %>%
  .[, .(ID, CHR, POS, A1, A2, A1_FREQ, OBS_CT, BETA, SE, T_STAT, P)]

fwrite(MTAGResults, file = paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWAS/output/MTAGinput/allChr_fpc",i,"_MTAG.txt"), sep="\t")

}) 


