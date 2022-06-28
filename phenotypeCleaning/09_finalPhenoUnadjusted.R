#!/usr/bin/env Rscript

## Run using R/4.1.2

library(data.table)
library(magrittr)
library(dplyr)
library(rlist)


imputedScans <- lapply(c(1:46), function(chunk) {
  paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/scansImputedWideFormat_chunk",chunk,".csv") %>%
  fread
  }) %>% rbindlist

pixels <- names(imputedScans)[!names(imputedScans) %in% c("patID", "eye", "visit", "scan", "refErr")]

## Function to generate phenotype
## If L and R available, take mean - otherwise take the selected single eye.

phenoAveraged <- imputedScans[, nScans := .N, by = "patID"] %>%
                 .[nScans == 2, eye := "both"]  %>%
                 .[, scan := NULL] %>%
                 .[, lapply(.SD, mean, na.rm=TRUE), .SDcols = c("refErr", pixels), by = c("patID", "eye", "visit")]  %>%
                 setnames(., "refErr", "meanRefErr")


## get other covariates - ie sex and age at OCT scan instance
file <- "/wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/app28541/phenoData/ukb41258.tab"
names <- fread(file, nrows = 0) %>%
 names

getIdx <- function(code) {
   names %like any% code %>% which
 }

cols <- getIdx(c("f.eid", "f.31.0.0", "f.21003.0%", "f.21003.1%",  "f.5270.%"))

covs <- fread(file, select = cols)

covsMerged <- phenoAveraged[, .(patID, visit)] %>%
 covs[., on = c(f.eid = "patID")] %>%
 .[, sex := f.31.0.0] %>%
 .[, age := case_when(visit == 0 ~ f.21003.0.0,
                     visit == 1 ~ f.21003.1.0)] %>%
 .[, device := case_when(visit == 0 ~ f.5270.0.0 %>% as.integer,
                     visit == 1 ~ f.5270.1.0)] %>%
 setnames(., "f.eid", "patID")

phenoOut <- covsMerged[, .(patID, visit, sex, age, device)] %>%
 .[phenoAveraged, on = c("patID", "visit")]

 fwrite(phenoOut, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/processedData/scansUnadjustedFinal.csv", sep = ",")
