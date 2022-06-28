#!/usr/bin/env Rscript

## Run using R/4.1.2

library(data.table)
library(magrittr)
library(dplyr)
library(parallel)
library(rlist)
library(ggplot2)
library(DescTools)


load("/stornext/Bioinf/data/lab_bahlo/projects/misc/retinalThickness/phenoExploratory/bestScansIndex20220504.RData")
refError <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/refractiveError.csv")
load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/gridMask.RData")

sampsIDs <- fread("/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/data/scanIDs.txt",
                col.names = c("id", "patientID"))

uniqueSamps <- sampsIDs[,id] %>% unique

coords <- missing0.1Final %>%
 reshape2::melt() %>%
 as.data.table %>%
 .[value == 0] %>%
 .[, value := NULL]


reshapeScans <- function(x) {

     sampIdx <- which(uniqueSamps==x)
     sampIdx %>% paste("Sample no.", .) %>% print

     scans <- bestScanIdx[[sampIdx]]

     sampScans <- lapply(scans, function(y) {

     # get ref error for eye
     eye <- ifelse(strsplit(y, "_")[[1]][2] == 21011, "L", "R")
     visit <- strsplit(y, "_")[[1]][3]

     refErr <- refError[id==x, get(paste0(eye,visit))]

     ## read in scan
    file <- paste0("/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/data/rawPerScan/",x,"/",y,"_scan.csv")

     if (file.exists(file))  {
     dt <- fread(file, col.names = c("patientID", "eye", "slice_index", "fovea_index", 0:255)) %>%
       .[, c("id", "laterality", "visit", "measure") := tstrsplit(patientID, "_")]

     fullScanMat <- dt[, c("slice_index", 0:255)] %>%
       as.matrix(., rownames="slice_index")

     # set -1 to NA
     fullScanMat[fullScanMat == -1] <- NA

     # reshape and merge with coords data.table
     trimmedScanLong <- fullScanMat %>%
      reshape2::melt() %>%
      as.data.table %>%
      .[coords, on = c("Var1", "Var2")] %>%
      .[, patID := x] %>%
      .[, eye := eye] %>%
      .[, visit := visit] %>%
      .[, scan := y] %>%
      .[, refErr := refErr]

      trimmedScanWide <-   trimmedScanLong  %>%
        dcast(., patID + eye + visit + scan + refErr ~ Var1 + Var2, value.var = "value")

      return(trimmedScanWide)
    }
  }) %>% rbindlist

  return(sampScans)
}

#
# scansDT <- lapply(uniqueSamps, reshapeScans) %>%
#  rbindlist
#
#
# fwrite(scansDT , "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/scansWideFormat.csv", sep = ",")

scansDT  <-fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/scansWideFormat.csv")



models <- lapply(pixels, function(y) {

  pixelIdx <- which(pixels==y)
  paste("pixel", pixelIdx) %>% print

  model <- lm(get(y) ~ refErr, data = scansDT, na.action=na.exclude)
  # return(model)

  beta <- model$coefficients["refErr"]
  t <- summary(model)$coefficients["refErr", "t value"]
  pVal <- summary(model)$coefficients["refErr", "Pr(>|t|)"]

  modelStats <- data.table(pixel =y,
    beta = beta,
    t = t,
    pVal = pVal)

    intercept <- model$coefficients["(Intercept)"]
    res <- residuals(model)

    out <- list(modelStats = modelStats,
      intercept = intercept,
      res = res)

    return(out)

})

names(models) <- pixels

save(models, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/reffErrorModels.RData")


pixelAssocsDT <- list.map(models, modelStats) %>%
  rbindlist %>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert=TRUE)]

adjustedPheno <- sapply(models, function(x) {

    adjustPheno <- x$intercept + x$res

    return(adjustPheno)

  }, USE.NAMES =T) %>%
  cbind(scansDT[, c("patID", "eye", "visit", "scan", "refErr")], .)



for(stat in c("beta", "t")) {

  plot <- ggplot(pixelAssocsDT , aes_string(x = "x", y = "y", fill = stat)) +
    geom_tile() +
    scale_fill_gradient2()

  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/plots/refErrorPixelwise_",stat,".png"))
  print(plot)
  dev.off()

}

minpval <- pixelAssocsDT[pVal>0, pVal] %>% min
pixelAssocsDTplot <- pixelAssocsDT %>%
  .[pVal==0, pVal := minpval] %>%
  .[, log10p := (-1)*log10(pVal)]


png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/plots/refErrorPixelwise_pVal.png")
ggplot(pixelAssocsDTplot , aes(x = x, y = y, fill = log10p)) +
  geom_tile()
dev.off()

fwrite(adjustedPheno , "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/scansRefErrAdjustedWideFormat.csv", sep = ",")
