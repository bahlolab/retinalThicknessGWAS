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

# scansDT <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/processedData/scansUnadjustedFinal.csv")
# pixels <-  names(scansDT)[!names(scansDT) %in% c("patID", "eye", "visit", "sex", "age", "device", "meanRefErr")]

# ids <- scansDT[,patID] 
# scansList <- lapply(ids, function(id) {
#   scansDT[patID==id, ..pixels] %>%
#    melt()  %>%
#   .[, c("y", "x") := tstrsplit(variable, "_", type.convert=TRUE)] %>%
#   dcast(., y~x, value.var = "value") %>%
#   as.matrix(., rownames = "y")
# })

# rm(scansDT)
# save(scansList, file =  "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/processedData/scansMatricesUnadjustedFinal.RData")

# ## summarise SD per pixel

# scansArray <- simplify2array(scansList)

# pixelSDs <- apply(scansArray, 1:2, sd, na.rm=T) %>%
# melt %>%
# as.data.table


# png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/pixelWiseSDs.png", width = 600, height = 600)
# ggplot(pixelSDs) +
#   geom_tile(aes(y = Var2, x = Var1, fill = value)) +
#   scale_fill_gradient2() +
#   scale_y_reverse() +
#   theme_bw() +
#   theme(legend.position = "bottom")
# dev.off()

# png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/pixelWiseVariance.png", width = 600, height = 600)
# ggplot(pixelSDs) +
#   geom_tile(aes(y = Var2, x = Var1, fill = value^2)) +
#   scale_fill_gradient2() +
#   scale_y_reverse() +
#   theme_bw() +
#   theme(legend.position = "bottom")
# dev.off()

# png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/pixelWiseVarianceHistograam.png", width = 600, height = 600)
# ggplot(pixelSDs, aes(x = value^2)) +
#   geom_histogram() +
#   theme_bw() +
#   theme(legend.position = "bottom")
# dev.off()

# pixelSDs <- pixelSDs[, Var50 := ifelse(value^2 > 50, 1, 0)] %>%
#  .[, Var55 := ifelse(value^2 > 55, 1, 0)] %>%
#  .[, Var60 := ifelse(value^2 > 60, 1, 0)] 

# png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/pixelWiseVarianceOver50.png", width = 600, height = 600)
# ggplot(pixelSDs) +
#   geom_tile(aes(y = Var2, x = Var1, fill = Var50)) +
#   scale_fill_gradient2() +
#   scale_y_reverse() +
#   theme_bw() +
#   theme(legend.position = "bottom")
# dev.off()

# png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/pixelWiseVarianceOver55.png", width = 600, height = 600)
# ggplot(pixelSDs) +
#   geom_tile(aes(y = Var2, x = Var1, fill = Var55)) +
#   scale_fill_gradient2() +
#   scale_y_reverse() +
#   theme_bw() +
#   theme(legend.position = "bottom")
# dev.off()

# png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/pixelWiseVarianceOver60.png", width = 600, height = 600)
# ggplot(pixelSDs) +
#   geom_tile(aes(y = Var2, x = Var1, fill = Var60)) +
#   scale_fill_gradient2() +
#   scale_y_reverse() +
#   theme_bw() +
#   theme(legend.position = "bottom")
# dev.off()


load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/processedData/scansMatricesUnadjustedFinal.RData")
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
             M = 50,
             # uniExpansions = list(list(type = "splines2D", k = c(10,10))),
             uniExpansions = list(list(type = "splines2Dpen",  k = c(12,12), parallel=T)),
             verbose = T)

save(pca, file="/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenoExploratory/working/cleanedScansFPCA_50fPCs_20221121.RData")


scans <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/processedData/scansUnadjustedFinal.csv")
pixels <-  names(scans)[!names(scans) %in% c("patID", "eye", "visit", "sex", "age", "device", "meanRefErr")]
coords <-tstrsplit(pixels, split = "_")
xMin <- coords[[2]] %>% min() %>% as.integer
yMin <- coords[[1]] %>% min() %>% as.integer


mean <- pca$meanFunction[[1]] %>%
  funData::as.data.frame(.) %>%
  as.data.table
  
png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcMeanFunction_50fPCs_20221121.png", width = 600, height = 600)
ggplot(mean) +
  geom_tile(aes(y = argvals1, x = argvals2, fill = X)) +
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

  return(plot)

})

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcFunctions_50fPCs_20221121.png", width = 1500, height = 1200)
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

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcScoreCorrelations_50fPCs_20221121.png", width = 1500, height = 1200)
reduce(fpcCorrPlots , `+`) %>%
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
