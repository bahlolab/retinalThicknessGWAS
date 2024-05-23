#!/usr/bin/env Rscript

## Run using R/4.1.2


## examine missingness patterns over the grid
## generate trimmed grids, with pixels (and rows/columns) with high missingness excluded

library(plotly)
library(data.table)
library(magrittr)
library(dplyr)
library(rlist)
library(parallel)
library(foreach)

mcf <- function(f) {
  function(...) {tryCatch({f(...)} , warning=function(w) { print(w) })
}}

getScans <- function(x) {

  which(uniqueScans==x) %>% paste("Scan no.", .) %>% print

  y <- sampsIDs[patientID == x, id]

  fullScanMat <- NULL

  # print(y)
  file <- paste0("/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/data/rawPerScan/",y,"/",x,"_scan.csv")

  if (file.exists(file))  {
    dt <- fread(file, col.names = c("patientID", "eye", "slice_index", "fovea_index", 0:255)) %>%
      .[, c("id", "laterality", "visit", "measure") := tstrsplit(patientID, "_")]

    fullScanMat <- dt[, c("slice_index", 0:255)] %>%
      as.matrix(., rownames="slice_index")


    # set -1 to NA
    fullScanMat[fullScanMat == -1] <- NA
  }
  return(fullScanMat)
}

plotSurface <- function(plotMat) {

  fig <- plot_ly(z = ~ plotMat) %>% add_surface(
    contours = list(
      z = list(
        show=TRUE,
        usecolormap=TRUE,
        highlightcolor="#ff0000",
        project=list(z=TRUE)
      )
    )
  )
  fig <- fig %>% layout(
    scene = list(
      camera=list(
        eye = list(x=1.87, y=0.88, z=-0.64)
      )
    )
  )
  return(fig)
}


sampsIDs <- fread("/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/data/scanIDs.txt",
                  col.names = c("id", "patientID"))

uniqueSamps <- sampsIDs[,id] %>% unique
fwrite(data.table(id = uniqueSamps), file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/uniqueSampsUnfiltered.txt", sep = "\t")

uniqueScans <- sampsIDs[, patientID] %>% unique
fwrite(data.table(id = uniqueSamps), file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/uniqueScansUnfiltered.txt", sep = "\t")

## if falls over when running in parallel, increase memory..
scansList <- mclapply(uniqueScans, getScans, mc.cores = 8)
names(scansList) <- uniqueScans

#scansList <- lapply(uniqueSamps, getScans)
save(scansList, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/scansListRaw.RData")
## Missingness pre pixel
# scansMissing <- mclapply(scansList, is.na, mc.cores = 6)
scansMissing <- lapply(scansList, is.na)

missingTotals <- Reduce('+', scansMissing)

## plot over grid
p <- plotSurface(missingTotals)
htmlwidgets::saveWidget(as_widget(p), "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/plots/gridMissingness.html")

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/plots/gridMissingness.png")
heatmap(missingTotals, Colv = NA, Rowv = NA, scale="column")
dev.off()


missingTotalsCapped <- missingTotals
missingTotalsCapped[missingTotals > 10000] <- 10000

p <- plotSurface(missingTotalsCapped)
htmlwidgets::saveWidget(as_widget(p), "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/plots/gridMissingnessCapped.html")

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/plots/gridMissingnessCapped.png")
heatmap(missingTotalsCapped, Colv = NA, Rowv = NA, scale="column")
dev.off()

## Exclude pixel with missingness great than n...
propMissing <- missingTotals / length(scansMissing)
missing0.05 <- ifelse(propMissing > 0.05, 1, 0)
missing0.1 <- ifelse(propMissing > 0.1, 1, 0)

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/plots/gridMissingnessProp0.05.png")
heatmap(missing0.05, Colv = NA, Rowv = NA, scale = "none")
dev.off()

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/plots/gridMissingnessProp0.10.png")
heatmap(missing0.1, Colv = NA, Rowv = NA, scale = "none")
dev.off()

## proportion of pixels removed, by row and by column, under the two thresholds
apply(missing0.05, 1, sum)
apply(missing0.1, 1, sum)
apply(missing0.05, 2, sum)
apply(missing0.1, 2, sum)


## remove whole rows / columns, if >50% pixels have >10% missing
removeRows <- apply(missing0.1, 1, sum) / ncol(missing0.1) > 0.5
removeCols <- apply(missing0.1, 2, sum) / nrow(missing0.1) > 0.5

missing0.1Final <- missing0.1
missing0.1Final[removeRows,] <- 1
missing0.1Final[,removeCols] <- 1

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/plots/gridMissingnessProp0.10Final.png")
heatmap(missing0.1Final, Colv = NA, Rowv = NA, scale = "none")
dev.off()

## Also remove pixels with high variance
# load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/scansListRaw.RData")
# load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/gridMask.RData")

scansArray <- simplify2array(scansList)
rm(scansList)

pixelSDs <- apply(scansArray, 1:2, sd, na.rm=T) %>%
melt %>%
as.data.table

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/plots/pixelWiseSDs.png", width = 600, height = 600)
ggplot(pixelSDs) +
  geom_tile(aes(y = Var2, x = Var1, fill = value)) +
  scale_fill_gradient2() +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "bottom")
dev.off()

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/plots/pixelWiseVariance.png", width = 600, height = 600)
ggplot(pixelSDs) +
  geom_tile(aes(y = Var2, x = Var1, fill = value^2)) +
  scale_fill_gradient2() +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "bottom")
dev.off()

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/plots/pixelWiseVarianceHistograam.png", width = 600, height = 600)
ggplot(pixelSDs, aes(x = value^2)) +
  geom_histogram() +
  theme_bw() +
  theme(legend.position = "bottom")
dev.off()

pixelSDs <- pixelSDs[, Var50 := ifelse(value^2 > 50, 1, 0)] %>%
 .[, Var55 := ifelse(value^2 > 55, 1, 0)] %>%
 .[, Var60 := ifelse(value^2 > 60, 1, 0)] 

fwrite(pixelSDs, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/pixelSDs.csv")

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/plots/pixelWiseVarianceOver50.png", width = 600, height = 600)
ggplot(pixelSDs) +
  geom_tile(aes(y = Var2, x = Var1, fill = Var50)) +
  scale_fill_gradient2() +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "bottom")
dev.off()

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/plots/pixelWiseVarianceOver55.png", width = 600, height = 600)
ggplot(pixelSDs) +
  geom_tile(aes(y = Var2, x = Var1, fill = Var55)) +
  scale_fill_gradient2() +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "bottom")
dev.off()

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/plots/pixelWiseVarianceOver60.png", width = 600, height = 600)
ggplot(pixelSDs) +
  geom_tile(aes(y = Var2, x = Var1, fill = Var60)) +
  scale_fill_gradient2() +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "bottom")
dev.off()

## Use threshold of over 55
highVar <- pixelSDs[, .(Var1, Var2, Var55)] %>%
    dcast(., Var1 ~ Var2) %>%
    as.matrix(.,  rownames = "Var1")

highVarMissing <- highVar + missing0.1Final

missing0.1Final <- ifelse(highVarMissing > 0, 1, 0)

save(missing0.1Final, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/gridMask.RData")

## plot high Var and missing0.1 Final on grid using ggplot and geom-tile

highVarPlot <- melt(highVar) %>%
  as.data.table %>%
  ggplot(.) +
  geom_tile(aes(y = Var1, x = Var2, fill = as.factor(value))) +
  # scale_fill_manual() +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "bottom") +
  ggtitle("Pixels with  variance > 55")

missingPlot <- melt(missing0.1Final) %>%
  as.data.table %>%
  ggplot(.) +
  geom_tile(aes(y = Var1, x = Var2, fill = as.factor(value))) +
  # scale_fill_manual() +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "bottom")  +
  ggtitle("Pixels excluded due to missingness")

library(patchwork)

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/plots/gridMissingnessProp0.10FinalHighVar.png", width = 1200, height = 600)
highVarPlot + missingPlot
dev.off()


## plot mean vs prop missing
pixelmeans <- apply(scansArray, 1:2, mean, na.rm=T) %>%
melt %>%
as.data.table %>%
setnames(., c("Var1", "Var2", "mean")) 
fwrite(pixelmeans, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/pixelMeans.csv")

load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/scansListRaw.RData")
scansMissing <- lapply(scansList, is.na)

missingTotals <- Reduce('+', scansMissing)

propMissing <- missingTotals / length(scansMissing)

pixelMissingProp <- propMissing %>%
melt %>%
as.data.table  %>%
setnames(., c("Var1", "Var2", "propMissing")) 

load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/gridMask.RData")
exclude <- missing0.1Final %>%
  melt %>%
  as.data.table %>%
  setnames(., c("Var1", "Var2", "exclude"))

missingnessDT <- pixelMissingProp[pixelmeans, on = c("Var1", "Var2")] %>%
  .[pixelSDs, on = c("Var1", "Var2")] %>%
  .[exclude, on = c("Var1", "Var2")] %>%
  .[exclude==0]

## plot propMissing vs mean, and propMissing vs SD (value) 
png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/plots/propMissingVsMean.png")
ggplot(missingnessDT) +
  geom_point(aes(x = mean, y = propMissing)) +
  theme_bw()
dev.off()

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/plots/propMissingVsSD.png")
ggplot(missingnessDT) +
  geom_point(aes(x = value, y = propMissing)) +
  theme_bw()
dev.off()



# trimmedScans <- lapply(scansList, function(x) { replace(x, missing0.1Final == 1, NA)})
# names(trimmedScans) <- uniqueScans
# save(trimmedScans, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/scansListTrimmed.RData")
