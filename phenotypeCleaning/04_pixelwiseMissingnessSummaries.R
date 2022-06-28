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
scansMissing <- mclapply(scansList, is.na, mc.cores = 6)

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

save(missing0.1Final, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/gridMask.RData")
# trimmedScans <- lapply(scansList, function(x) { replace(x, missing0.1Final == 1, NA)})
# names(trimmedScans) <- uniqueScans
# save(trimmedScans, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/scansListTrimmed.RData")
