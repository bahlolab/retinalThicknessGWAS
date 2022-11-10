#!/usr/bin/env Rscript

library(data.table)
library(magrittr)
library(tidyverse)


covars <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/covariates_doubleIDs.txt")    
fpcs <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWAS/output/FPCphenotypes_doubleIDs.txt")

 newCovars <- covars[fpcs[, .(FID, IID, fpc1)], on = c("FID", "IID")]

 fwrite(newCovars, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/output/covariates_withfpc1_doubleIDs.txt", sep = "\t", na = "NA", quote = F)
 