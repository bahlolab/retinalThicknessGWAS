## Run using R/4.1.2


library(data.table)
library(magrittr)
library(dplyr)
library(here)

hits <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/cojoAnalysis/output/bonfSigLociWithSecondary.csv") %>%
.[, c("slice", "x") := tstrsplit(pixel, "_", type.convert=TRUE)] 

hitsPrimary <- hits[, .(locusID, CHR, ID, pixel, slice)] %>%
    unique 

hitsSecondary <-  hits[, .(locusID, CHR, SNP_secondary, pixel, slice)] %>%
.[SNP_secondary!=""] %>%
setnames(., "SNP_secondary", "ID") 


hitsOut <- rbind(hitsPrimary, hitsSecondary) %>%
    setkey(., "locusID") 
    
  
fwrite(hitsOut,  file = here("loci.txt"), sep = "\t",  quote = F, col.names = F)