#!/usr/bin/env Rscript

## Run using R/4.1.2

library(data.table)
library(magrittr)
library(dplyr)
library(rlist)
library(DescTools)

bestScansDT <- fread("/stornext/Bioinf/data/lab_bahlo/projects/misc/retinalThickness/phenoExploratory/bestScansIDT20221118.txt")

imputedScans <- lapply(c(1:46), function(chunk) {
  paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/scansImputedWideFormat_chunk",chunk,".csv") %>%
  fread %>%
    .[scan %in% bestScansDT[,scan]]
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

## Final filtering of extreme values / outliers

pixels <-  names(phenoOut)[!names(scansDT) %in% c("patID", "eye", "visit", "sex", "age", "device", "meanRefErr")]

minimums <- apply(phenoOut[, ..pixels], 1, min, na.rm = T)
maximums <- apply(phenoOut[, ..pixels], 1, max, na.rm = T)
means <- apply(phenoOut[, ..pixels], 1, mean, na.rm = T)
sds <- apply(phenoOut[, ..pixels], 1, sd, na.rm = T)

summaries <- data.table(patID = phenoOut[,patID],
                        minimum = minimums,
                        maximum = maximums,
                        mean = means,
                        sd = sds) 
summaries %>% summary
summaries <- summaries[, outlier := case_when( minimum < 30 ~ 1,
                            maximum > 150 ~ 1,
                            sd > 15 ~ 1,
                            T ~ 0)]
remove <- summaries[outlier==1, patID]

phenoOut <- phenoOut[!patID %in% remove]

