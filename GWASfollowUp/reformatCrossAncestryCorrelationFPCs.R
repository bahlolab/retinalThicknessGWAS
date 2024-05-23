#!/usr/bin/env Rscript

library(data.table)
library(magrittr)
library(here)
library(dplyr)

for(ANC in c("EUR", "CSA", "AFR")) {

  for(fpc in c(1:6)) {

  outData <- lapply(c(1:22, "X"), function(chr) {


  if(ANC == "EUR") {

    file <- here::here("fpcResults", paste0("chr",chr,ANC,".fpc",fpc,".glm.linear"))
    
    print(paste0("Processing ", file))
    
    data <- fread(file) %>%
    setnames(., c("CHR", "BP", "SNP", "REF", "ALT", "A1", "A1_CT", "freq", "TEST", "N", "b", "se", "Tstat", "p", "Err")) %>%
    .[, A2 := case_when(A1 == REF ~ ALT, 
                        A1 == ALT ~ REF)] %>%
    .[freq %between% c(0.05, 0.95)] %>%
    .[, .(CHR, BP, SNP, A1, A2,  N,  b, se, p)] %>%
    setnames(., c("chr", "pos", "SNP", "A1", "A2", "N", "beta", "SE", "p-value"))

    if(chr == 6) {
      data <- data %>%
      .[ !(pos %between% c(28477797, 33448354))] 
    }

  }  else {
  file <-  here::here("fpcResults", paste0("chr",chr,"_",ANC,".fpc",fpc,".glm.linear"))

    print(paste0("Processing ", file))
    data <- fread(file) %>%
    setnames(., c("CHR", "BP", "SNP", "REF", "ALT", "PROVISIONAL_REF", "A1", "OMITTED",  "A1_CT", "freq", "TEST", "N", "b", "se", "Tstat", "p", "Err")) %>%
    .[, A2 := case_when(A1 == REF ~ ALT, 
                        A1 == ALT ~ REF)] %>%
    .[freq %between% c(0.05, 0.95)] %>%
    .[, .(CHR, BP, SNP, A1, A2,  N,  b, se, p)] %>%
    setnames(., c("chr", "pos", "SNP", "A1", "A2", "N", "beta", "SE", "p-value"))

        if(chr == 6) {
      data <- data %>%
      .[ !(pos %between% c(28477797, 33448354))] 
    }

  }
    return(data)

    }) %>%
    rbindlist

      outName <- paste0("allChr.fpc",fpc,"_",ANC,"popcornFormat.txt") %>%
  here::here("popcornInFiles", .)
    
  print(paste0("Writing to ", outName))
  fwrite(outData, outName, sep = "\t")


  }
}

