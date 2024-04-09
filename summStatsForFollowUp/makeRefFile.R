## Run using R/4.1.2


library(data.table)
library(magrittr)
library(dplyr)
library(here)

hits <- fread("/stornext/Bioinf/data/lab_bahlo/projects/misc/retinalThickness/GWAS/finalResultsEUR/gwSigLociSummary.csv") %>%
.[, c("slice", "x") := tstrsplit(pixel, "_", type.convert=TRUE)] %>%
.[BonferroniSig == "Y"] %>%
.[, .(locusID, CHR, ID, pixel, slice)]

fwrite(hits,  file = here("loci.txt"), sep = "\t",  quote = F, col.names = F)