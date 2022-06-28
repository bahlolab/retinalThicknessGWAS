#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## Run using R/4.1.2

library(data.table)
library(magrittr)
library(dplyr)
# library(parallel)
# library(doParallel)
# library(iterators)
# library(foreach)
library(rlist)
library(ggplot2)

chunk <- args[1] %>% as.numeric

print(paste("chunk",chunk))

# chunkedPheno <- iter(pheno, by="row", chunksize=1000)

start <- ((chunk - 1) * 2000) + 1
end <- chunk * 2000

   # nCores <- 8
   #
   # cluster <- makeCluster(nCores)
   # registerDoParallel(cluster)

 pheno <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/scansWideFormat.csv")
 sampsIDs <- fread("/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/data/scanIDs.txt",
                 col.names = c("id", "patientID"))

 pixels <- names(pheno)[!names(pheno) %in% c("patID", "eye", "visit", "scan", "refErr")]
 uniqueScans <- pheno[,scan] %>% unique


set.seed(3467)

chunkedPheno <- pheno[start : end] %>%
  .[!is.na(patID)]

rm(pheno)

imputeScans <- function(scanidx) {

     scanIdx <- which(uniqueScans==scanidx)
     scanIdx %>% paste("Scan no.", .) %>% print

     scanDT <- chunkedPheno[scan == scanidx, ..pixels] %>%
       reshape2::melt() %>%
       as.data.table %>%
       .[, c("x", "y") := tstrsplit(variable, "_", type.convert=TRUE)]

      if(nrow(scanDT[is.na(value)]) > 0) {
        mod <- mgcv::bam(value ~ te(x, y, bs =  "ps",  k = c(12,12)), data = scanDT, method = "REML", na.action = "na.exclude")

        imputedDF <- scanDT %>%
          .[, pred := predict(mod, scanDT, type = "response")] %>%
          # .[, resdiual := residuals.gam(mod)] %>%
          .[, imputed := ifelse(is.na(value), pred, value)]

        imputed <- imputedDF[, imputed]

      } else {

        imputed <- scanDT[,  value]

      }

    names(imputed) <- pixels

    imputedScan <- setDT(as.list(imputed)) %>%
     cbind(chunkedPheno[scan == scanidx , c("patID", "eye", "visit", "scan", "refErr")],  .)

   return(imputedScan)

  }



chunkedUniqueScans <- chunkedPheno[,scan] %>% unique

imputedScans <- lapply(chunkedUniqueScans, imputeScans) %>%
rbindlist



fwrite(imputedScans, file = paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/scansImputedWideFormat_chunk",chunk,".csv"), sep = ",")














## old parallelisation...
## doesn't work

# imputedScans <- foreach(scan = chunkedPheno,
#                         .combine = rbind,
#                         .packages = c("magrittr", "data.table", "mgcv", "dplyr")) %dopar% {
#
#   # scanIdx <- which(uniqueScans==scan$scan)
#   # scanIdx %>% paste("Scan no.", .) %>% print
#
#    scanDT <- scan[, ..pixels] %>%
#      reshape2::melt() %>%
#      as.data.table %>%
#      .[, c("x", "y") := tstrsplit(variable, "_", type.convert=TRUE)]
#
#     if(nrow(scanDT[is.na(value)]) > 0) {
#       mod <- mgcv::bam(value ~ te(x, y, bs =  "ps",  k = c(12,12)), data = scanDT, method = "REML", na.action = "na.exclude")
#
#       imputedDF <- scanDT %>%
#         .[, pred := predict(mod, scanDT, type = "response")] %>%
#         # .[, resdiual := residuals.gam(mod)] %>%
#         .[, imputed := ifelse(is.na(value), pred, value)]
#
#       imputed <- imputedDF[, imputed]
#
#     } else {
#
#       imputed <- scanDT[,  value]
#
#     }
#
#   names(imputed) <- pixels
#
#   imputedScan <- setDT(as.list(imputed)) %>%
#    cbind(scan[ , c("patID", "eye", "visit", "scan", "refErr")],  .)
#
#  return(imputedScan)
#
# }
#
#  end_time <- Sys.time()
#  print(end_time - start_time)
#  timeV1 <- end_time - start_time
#

# stopCluster(cluster)
