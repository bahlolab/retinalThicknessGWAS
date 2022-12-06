#!/usr/bin/env Rscript

library(MFPCA)
library(data.table)
library(magrittr)
library(dplyr)
library(rlist)
library(doParallel)
library(foreach)
library(abind)
library(ggplot2)
library(patchwork)
library(purrr)

nCores <- 6
cluster <- makeCluster(nCores)
doParallel::registerDoParallel(cluster)

# detectedCores <- detectCores()
#
# print(paste("Sense check - there are",detectedCores," cores!"))
set.seed(3467)

load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/scansMatricesUnadjustedFinalfPCexclusions.Data")
scansArray <- simplify2array(scansList) %>%
  aperm(., c(3, 2, 1))

# rows <- c(2:113) %>% as.character()
# cols <- c(9:250) %>% as.character()

# # dim2idx <- dimnames(scansArray)[2] %in% rows %>% which
# # dim3idx <- dimnames(scansArray)[3] %in% cols %>% which

# scansArrayFilt <- asub(scansArray, idx = list(rows, cols), dims = c(2,3))

y <- dim(scansArray)[2]
z <- dim(scansArray)[3]

domain <- list(c(1:y), c(1:z))

scansFunData <- funData(domain, scansArray)

scansMultiFunData <- multiFunData(list(scansFunData))

## tidy up before running pca
rm(scansList, scansArray, scansFunData)

pca <- MFPCA(scansMultiFunData,
             M = 100,
             # uniExpansions = list(list(type = "splines2D", k = c(10,10))),
             uniExpansions = list(list(type = "splines2Dpen",  k = c(12,12), parallel=T)),
             verbose = T)

save(pca, file="/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/cleanedScansFPCAround2_100fPCs_20221124.RData")


scans <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/scansUnadjustedFinalfPCexclusions.csv")
pixels <-  names(scans)[!names(scans) %in% c("patID", "eye", "visit", "sex", "age", "device", "meanRefErr")]
coords <-tstrsplit(pixels, split = "_")
xMin <- coords[[2]] %>% min() %>% as.integer
yMin <- coords[[1]] %>% min() %>% as.integer

mean <- pca$meanFunction[[1]] %>%
  funData::as.data.frame(.) %>%
  as.data.table
  
png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcMeanFunction_round2_100fPCs_20221124.png", width = 600, height = 600)
ggplot(mean) +
  geom_tile(aes(x = argvals1, y = argvals2, fill = X)) +
  scale_fill_gradient2() +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "bottom")
dev.off()


mask <- mean[,X] %>% is.na %>% which

fpcPlots <- lapply(c(1:20), function(i) {

   fpc <- pca$functions[[1]] %>%
    funData::as.data.frame(.) %>%
    as.data.table %>%
    .[obs==i] %>%
    .[!mask] %>%
    .[, x := argvals1 + (xMin - 1)] %>%
    .[, y := argvals2 + (yMin - 1)] 

  plot <- ggplot(fpc) +
    geom_tile(aes(x = x, y = y, fill = X)) +
    scale_fill_gradient2() +
    scale_y_reverse() +
    theme_bw() +
    theme(legend.position = "none")
  # theme(legend.position = "bottom")

  return(plot)

})

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcFunctions_round2_100fPCs_20221124.png", width = 1500, height = 1200)
reduce(fpcPlots , `+`) %>%
  print
dev.off()


reorderIdx <- match(rownames(pca$scores), as.character(scans[, patID]))
scans <- scans[reorderIdx]

fpcCorrPlots <- lapply(c(1:20), function(i) {

fpcCorr <- cor(scans[, ..pixels], pca$scores[,i]) %>%
  as.data.table(keep.rownames = T) %>%
  .[, c("y", "x") := tstrsplit(rn, "_", type.convert=TRUE)]

setnames(fpcCorr, "V1", "cor")

plot <- ggplot(fpcCorr) +
  geom_tile(aes(x = x, y = y, fill = cor)) +
  scale_fill_gradient2() +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "none")

return(plot)
})

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcScoreCorrelations_round2_100fPCs_20221124.png", width = 1500, height = 1200)
reduce(fpcCorrPlots , `+`) %>%
print
dev.off()



load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenoExploratory/working/cleanedScansFPCA_100fPCs_20221121.RData")
pca1 <- pca
scoresDT1 <- pca1$scores  %>%
as.data.table(keep.rownames = T) %>%
setnames(., c("patID", paste0("v1FPC",1:100))) 
load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/cleanedScansFPCAround2_100fPCs_20221124.RData")
pca2 <- pca
scoresDT2 <- pca2$scores  %>%
as.data.table(keep.rownames = T) %>%
setnames(., c("patID", paste0("v2FPC",1:100)))

scoresDT <- scoresDT1[scoresDT2, on  = "patID"]

verCorrPlots <- lapply(c(1:25), function(x) {
  
 ggplot(scoresDT, aes_string(paste0("v1FPC",x), paste0("v2FPC",x))) +
    geom_point() 
  
  })

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcScores_round1vs2_100fPCs.png", width = 1500, height = 1500)
reduce(verCorrPlots , `+`) %>%
  print
dev.off()


# scoresDT <- pca$scores  %>%
#    as.data.table(keep.rownames = T) %>%
#    setnames(., c("patID", paste0("fpc",1:100)))


# exclude <- lapply(c(1:100), function(i){

#   col <- paste0("fpc",i)
#   vals <- scoresDT[, get(col)] %>% as.vector

#   m <- vals  %>% mean
#   sd <- vals  %>% sd

#   exclude <- ifelse(vals %between% c((m-(10*sd)), (m+(10*sd))), 0, 1)

# })

# excludeAll <- list.cbind(exclude) %>%
#   as.data.table %>%
#   cbind(scoresDT[, "patID"], .) %>%
#   .[, excSum := rowSums(.SD, na.rm=T), .SDcols = paste0("V",c(1:100))] %>%
#   .[, exclude := ifelse(excSum >0, 1, 0)]

# exclIDs <- excludeAll[exclude==1, patID]

# scoresDT <- scoresDT[!patID %in% exclIDs]

# fwrite(scoresDT, file = "processedData/fpcPhenotypes.txt")


parallel::stopCluster(cluster)
