#!/usr/bin/env Rscript

library(MFPCA)
library(data.table)
library(magrittr)
library(dplyr)
library(rlist)
library(doParallel)
library(foreach)
library(abind)

nCores <- 6
cluster <- makeCluster(nCores)
doParallel::registerDoParallel(cluster)

# detectedCores <- detectCores()
#
# print(paste("Sense check - there are",detectedCores," cores!"))
set.seed(3467)


load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenoExploratory/working/cleanFinalScans.RData")

scansArray <- simplify2array(scansList) %>%
  aperm(., c(3, 2, 1))

rows <- c(2:113) %>% as.character()
cols <- c(9:250) %>% as.character()

# dim2idx <- dimnames(scansArray)[2] %in% rows %>% which
# dim3idx <- dimnames(scansArray)[3] %in% cols %>% which

scansArrayFilt <- asub(scansArray, idx = list(rows, cols), dims = c(2,3))

y <- dim(scansArrayFilt)[2]
z <- dim(scansArrayFilt)[3]

domain <- list(c(1:y), c(1:z))

scansFunData <- funData(domain, scansArrayFilt)

scansMultiFunData <- multiFunData(list(scansFunData))

## tidy up before running pca
rm(scansList, scansArray, scansArrayFilt, scansFunData)

pca <- MFPCA(scansMultiFunData,
             M = 100,
             # uniExpansions = list(list(type = "splines2D", k = c(10,10))),
             uniExpansions = list(list(type = "splines2Dpen",  k = c(12,12), parallel=T)),
             verbose = T)

save(pca, file="/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenoExploratory/working/cleanedScansFPCA_20220927.RData")

scoresDT <- pca$scores  %>%
   as.data.table(keep.rownames = T) %>%
   setnames(., c("patID", paste0("fpc",1:100)))


exclude <- lapply(c(1:100), function(i){

  col <- paste0("fpc",i)
  vals <- scoresDT[, get(col)] %>% as.vector

  m <- vals  %>% mean
  sd <- vals  %>% sd

  exclude <- ifelse(vals %between% c((m-(10*sd)), (m+(10*sd))), 0, 1)

})

excludeAll <- list.cbind(exclude) %>%
  as.data.table %>%
  cbind(scoresDT[, "patID"], .) %>%
  .[, excSum := rowSums(.SD, na.rm=T), .SDcols = paste0("V",c(1:100))] %>%
  .[, exclude := ifelse(excSum >0, 1, 0)]

exclIDs <- excludeAll[exclude==1, patID]

scoresDT <- scoresDT[!patID %in% exclIDs]

fwrite(scoresDT, file = "processedData/fpcPhenotypes.txt")


parallel::stopCluster(cluster)
