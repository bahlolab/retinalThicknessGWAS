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


load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenoExploratory/working/cleanFinalScans.RData")

## summarise SD per pixel
scansArray <- simplify2array(scansList)

pixelSDs <- apply(scansArray, 1:2, sd, na.rm=T) %>%
melt %>%
as.data.table

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/pixelWiseSDs.png", width = 600, height = 600)
ggplot(pixelSDs) +
  geom_tile(aes(y = Var2, x = Var1, fill = value)) +
  scale_fill_gradient2() +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "bottom")
dev.off()

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/pixelWiseVariance.png", width = 600, height = 600)
ggplot(pixelSDs) +
  geom_tile(aes(y = Var2, x = Var1, fill = value^2)) +
  scale_fill_gradient2() +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "bottom")
dev.off()

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/pixelWiseVarianceHistograam.png", width = 600, height = 600)
ggplot(pixelSDs, aes(x = value^2)) +
  geom_histogram() +
  theme_bw() +
  theme(legend.position = "bottom")
dev.off()

pixelSDs <- pixelSDs[, Var50 := ifelse(value^2 > 50, 1, 0)] %>%
 .[, Var55 := ifelse(value^2 > 55, 1, 0)] %>%
 .[, Var60 := ifelse(value^2 > 60, 1, 0)] 

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/pixelWiseVarianceOver50.png", width = 600, height = 600)
ggplot(pixelSDs) +
  geom_tile(aes(y = Var2, x = Var1, fill = Var50)) +
  scale_fill_gradient2() +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "bottom")
dev.off()

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/pixelWiseVarianceOver55.png", width = 600, height = 600)
ggplot(pixelSDs) +
  geom_tile(aes(y = Var2, x = Var1, fill = Var55)) +
  scale_fill_gradient2() +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "bottom")
dev.off()

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/pixelWiseVarianceOver60.png", width = 600, height = 600)
ggplot(pixelSDs) +
  geom_tile(aes(y = Var2, x = Var1, fill = Var60)) +
  scale_fill_gradient2() +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "bottom")
dev.off()

## Use threshold of over 50
highVar <- pixelSDs[, .(Var1, Var2, Var50)] %>%
    dcast(., Var1 ~ Var2) %>%
    as.matrix(.,  rownames = "Var1")

## filter scans
scansListFilt <- lapply(scansList, function(scan) {

scan <- replace(scan, highVar == 1, NA) 
class(scan) <- "numeric"
return(scan)

})


scansArrayFilt <- simplify2array(scansListFilt) %>%
  aperm(., c(3, 2, 1))


## rough manual row/col removal test
# rows <- c(2:113) %>% as.character()
# cols <- c(9:250) %>% as.character()

# # dim2idx <- dimnames(scansArray)[2] %in% rows %>% which
# # dim3idx <- dimnames(scansArray)[3] %in% cols %>% which

# scansArrayFilt <- asub(scansArray, idx = list(rows, cols), dims = c(2,3))

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

save(pca, file="/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenoExploratory/working/cleanedScansFPCA_var50Removed_20221110.RData")




scans <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/processedData/scansUnadjustedFinal.csv")
pixels <-  names(scans)[!names(scans) %in% c("patID", "eye", "visit", "sex", "age", "device", "meanRefErr")]

eigen <- data.table(fpc = c(1:100), val=pca$values)
# ggplot(eigen, aes(x = fpc, y=val)) +
#   geom_line() +
#   geom_point()

mean <- pca$meanFunction[[1]] %>%
  funData::as.data.frame(.) %>%
  as.data.table
  
png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcMeanFunction_var50Removed.png", width = 600, height = 600)
ggplot(mean) +
  geom_tile(aes(y = argvals1, x = argvals2, fill = X)) +
  scale_fill_gradient2() +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "bottom")
dev.off()



fpcPlots <- lapply(c(1:20), function(i) {

   fpc <- pca$functions[[1]] %>%
    funData::as.data.frame(.) %>%
    as.data.table %>%
    .[obs==i]

  plot <- ggplot(fpc) +
    geom_tile(aes(y = argvals1, x = argvals2, fill = X)) +
    scale_fill_gradient2() +
    scale_y_reverse() +
    theme_bw() +
    theme(legend.position = "none")
  # theme(legend.position = "bottom")

  return(plot)

})

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcFunctions_var50Removed.png", width = 1500, height = 1200)
reduce(fpcPlots , `+`) %>%
  print
dev.off()



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

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcScoreCorrelations_var50Removed.png", width = 1500, height = 1200)
reduce(fpcCorrPlots , `+`) %>%
print
dev.off()

