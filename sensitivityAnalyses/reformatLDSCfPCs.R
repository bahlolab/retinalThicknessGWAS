#!/usr/bin/env Rscript

library(data.table)
library(magrittr)
library(here)

refs <- lapply(c(1:22), function(chr) {
  here(paste0("chr",chr,"SNPinfo.txt")) %>%
  fread %>%
  setnames(., c("SNP", "CHR", "BP", "A1", "A2", "freq", "GENE"))
})

lapply(c(1:6), function(fpc) {

outName <- paste0("fPC",fpc, "_ldscFormat.txt") %>%
here::here("ldscInFiles", "fPCs",  .)

fpcResult <- lapply(c(1:22), function(chr) {

ref <- refs[[chr]]

data <- here::here("fpcResults",  paste0("chr",chr,"EUR.fpc",fpc,".glm.linear")) %>%
  fread(.) %>%
  setnames(., c("CHR", "BP", "SNP", "REF", "ALT", "A1", "A1_CT", "A1_FREQ", "TEST", "OBS_CT", "b", "se", "Tstat", "p", "ERRCODE")) %>%
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
fwrite(fpcResult, outName, sep = " ")

rm(fpcResult)

## smoking
outName <- paste0("fPC",fpc, "_smoking_ldscFormat.txt") %>%
here::here("ldscInFiles", "fPCs",  .)

fpcResult <- lapply(c(1:22), function(chr) {

ref <- refs[[chr]]

data <- here::here("results",  paste0("chr",chr,"_smoking.fpc",fpc,".glm.linear")) %>%
  fread(.) %>%
  setnames(., c("BP", "SNP", "REF", "ALT", "PROVISIONAL_REF", "A1", "OMITTED", "A1_FREQ",  "b", "se", "Tstat", "p")) %>%
  .[ref, on = c("SNP", "A1")] %>%
  .[, N := 42971] %>%
  .[, .(SNP, A1, A2, N, freq, p, b)] %>%
  setnames(., "freq", "FRQ") %>%
  setnames(., "b", "BETA") %>%
  setnames(., "p", "P")

  return(data)

}) %>%
rbindlist

print(paste0("Writing to ", outName))
fwrite(fpcResult, outName, sep = " ")

rm(fpcResult)

## smoking
outName <- paste0("fPC",fpc, "_noSurgery_ldscFormat.txt") %>%
here::here("ldscInFiles", "fPCs",  .)

fpcResult <- lapply(c(1:22), function(chr) {

ref <- refs[[chr]]

data <- here::here("results",  paste0("chr",chr,"_noSurgery.fpc",fpc,".glm.linear")) %>%
  fread(.) %>%
  setnames(., c("BP", "SNP", "REF", "ALT", "PROVISIONAL_REF", "A1", "OMITTED", "A1_FREQ",  "b", "se", "Tstat", "p")) %>%
  .[ref, on = c("SNP", "A1")] %>%
  .[, N := 43098] %>%
  .[, .(SNP, A1, A2, N, freq, p, b)] %>%
  setnames(., "freq", "FRQ") %>%
  setnames(., "b", "BETA") %>%
  setnames(., "p", "P")

  return(data)

}) %>%
rbindlist

print(paste0("Writing to ", outName))
fwrite(fpcResult, outName, sep = " ")

rm(fpcResult)


})




