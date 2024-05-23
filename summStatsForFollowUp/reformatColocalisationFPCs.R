#!/usr/bin/env Rscript

library(data.table)
library(magrittr)
library(here)
library(dplyr)

for(fpc in c(1:6)) {

  for(chr in c(1:22, "X")) {

file <-  here::here("fpcResults", paste0("chr",chr,"EUR.fpc",fpc,".glm.linear"))

outName <- paste0("chr",chr,".fpc",fpc,"_colocFormat.txt") %>%
here::here("colocInFiles", .)


  print(paste0("Processing ", file))
  data <- fread(file) %>%
  setnames(., c("CHR", "BP", "SNP", "REF", "ALT", "A1",  "A1_CT", "freq", "TEST", "N", "b", "se", "Tstat", "p", "Err")) %>%
  .[, A2 := case_when(A1 == REF ~ ALT, 
                      A1 == ALT ~ REF)] %>%
  .[, .(CHR, BP, SNP, A1, A2, freq, b, se, p, N)]

  print(paste0("Writing to ", outName))
  fwrite(data, outName, sep = " ")

  }
}


