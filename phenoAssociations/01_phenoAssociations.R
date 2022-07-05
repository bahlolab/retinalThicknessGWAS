#!/usr/bin/env Rscript

## Run using R/4.1.2

library(data.table)
library(magrittr)
library(dplyr)
library(rlist)
library(ggplot2)
library(DescTools)


dataDir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenoAssociations/rawData/"
outDir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenoAssociations/output/"
scansDT  <-fread(paste0(dataDir,"scansUnadjustedFinal.csv"))
load(paste0(dataDir,"Data_for_retinal_areas.RData"))

pixels <- names(scansDT)[!names(scansDT) %in% c("patID", "visit", "eye", "sex", "age", "device", "meanRefErr")]


# linkage <- fread(paste0(dataDir,"idLinkage.txt"))
# scansDTlinked <- linkage[scansDT, on = c("patIDhda" = "patID")]

file <- paste0(dataDir,"ukb41258.tab")
names <- fread(file, nrows = 0) %>%
  names

getIdx <- function(code) {
    names %like any% code %>% which
  }

cols <- getIdx(c("f.eid", "f.48.0.0", "f.48.1.0", "f.49.0.0", "f.49.1.0",
"f.50.0.0", "f.50.1.0", "f.51.0.0", "f.51.1.0", "f.93.0.%", "f.93.1.%",
"f.94.0.%", "f.94.1.%", "f.4079.0.%", "f.4079.1.%", "f.4080.0.%", "f.4080.1.%",
 "f.6177.0.%", "f.6177.1.%", "f.21002.%", "f.21001.%"))

covs <- fread(file, select = cols)

covsMerged <- scansDT[, .(patID, visit)] %>%
  covs[., on = c(f.eid = "patID")] %>%
  .[, waist := case_when(visit == 0 ~ f.48.0.0,
                      visit == 1 ~ f.48.1.0)] %>%
  .[, hip := case_when(visit == 0 ~ f.49.0.0,
                      visit == 1 ~ f.49.1.0)] %>%
  .[, standHeight := case_when(visit == 0 ~ f.50.0.0,
                      visit == 1 ~ f.50.1.0)] %>%
  .[, sitHeight := case_when(visit == 0 ~ f.51.0.0,
                      visit == 1 ~ f.51.1.0)] %>%
  .[, weight := case_when(visit == 0 ~ f.21002.0.0,
                      visit == 1 ~ f.21002.1.0)] %>%
  .[, bmi := case_when(visit == 0 ~ f.21001.0.0,
                      visit == 1 ~ f.21001.1.0)] %>%
  .[, SBP0 := rowMeans(.SD, na.rm=T), .SDcols = c("f.93.0.0", "f.93.0.1", "f.4080.0.0", "f.4080.0.1"), by =  f.eid]  %>%
  .[, SBP1 := rowMeans(.SD, na.rm=T), .SDcols = c("f.93.1.0", "f.93.1.1", "f.4080.1.0", "f.4080.1.1"), by =  f.eid]  %>%
  .[, DBP0 := rowMeans(.SD, na.rm=T), .SDcols = c("f.94.0.0", "f.94.0.1", "f.4079.0.0", "f.4079.0.1"), by =  f.eid]  %>%
  .[, DBP1 := rowMeans(.SD, na.rm=T), .SDcols = c("f.94.1.0", "f.94.1.1", "f.4079.1.0", "f.4079.1.1"), by =  f.eid]  %>%
  .[, rawSBP := case_when(visit == 0 ~ SBP0,
                      visit == 1 ~ SBP1)] %>%
  .[, rawDBP := case_when(visit == 0 ~ DBP0,
                      visit == 1 ~ DBP1)] %>%
  .[, BPmedication := case_when(visit == 0 & (f.6177.0.0 %in% 2 | f.6177.0.1 %in% 2 | f.6177.0.2 %in% 2) ~ 1,
                      visit == 0 & !(f.6177.0.0 %in% 2 | f.6177.0.1 %in% 2 | f.6177.0.2 %in% 2) ~ 0,
                      visit == 1 & (f.6177.0.0 %in% 2 | f.6177.0.1 %in% 2 | f.6177.0.2 %in% 2) ~ 1,
                      visit == 1 & !(f.6177.0.0 %in% 2 | f.6177.0.1 %in% 2 | f.6177.0.2 %in% 2) ~ 0)] %>%
  .[, adjustedSBP := case_when(visit == 0 & BPmedication == 1 ~ SBP0 + 15,
                      visit == 0 & BPmedication == 0 ~ SBP0,
                      visit == 1 & BPmedication == 1 ~ SBP1 + 15,
                      visit == 1 & BPmedication ==0 ~ SBP1)] %>%
  .[, adjustedDBP := case_when(visit == 0 & BPmedication ==1 ~ DBP0 + 10,
                      visit == 0 & BPmedication == 0 ~ DBP0,
                      visit == 1 & BPmedication == 1 ~ DBP1 + 10,
                      visit == 1 & BPmedication ==0 ~ DBP1)] %>%
  setnames(., "f.eid", "patID")

phenoOut <- covsMerged[, .(patID, visit, waist, hip, standHeight, sitHeight, weight, bmi, rawSBP, rawDBP, adjustedSBP, adjustedDBP, BPmedication)] %>%
  .[scansDT, on = c("patID", "visit")] %>%
  .[, ageSq := age^2]

phenos <- c("sex", "age", "ageSq", "waist", "hip", "standHeight", "sitHeight", "weight", "bmi", "rawSBP", "rawDBP", "adjustedSBP", "adjustedDBP", "BPmedication")

models <- lapply(pixels, function(y) {

  pixelIdx <- which(pixels==y)
  paste("pixel", pixelIdx) %>% print

  modelStats <- lapply(phenos, function(x) {
if(x == "ageSq") {
  model <- lm(get(y) ~ age + ageSq, data = phenoOut, na.action=na.exclude)

  beta <- model$coefficients["ageSq"]
  t <- summary(model)$coefficients["ageSq", "t value"]
  pVal <- summary(model)$coefficients["ageSq", "Pr(>|t|)"]

} else {
  model <- lm(get(y) ~ get(x), data = phenoOut, na.action=na.exclude)

  beta <- model$coefficients["get(x)"]
  t <- summary(model)$coefficients["get(x)", "t value"]
  pVal <- summary(model)$coefficients["get(x)", "Pr(>|t|)"]
}

  out<- data.table(pixel = y,
    beta = beta,
    t = t,
    pVal = pVal)

    return(out)

  })
  names(modelStats) <- phenos
  return(modelStats)
})
names(models) <- pixels

save(models, file = paste0(outDir,"phenoModels.RData"))

# load(paste0(outDir,"phenoModels.RData"))

lapply(phenos, function(trait) {

pixelAssocsDT <- list.map(models, get(trait)) %>%
  rbindlist %>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert=TRUE)]


for(stat in c("beta", "t")) {

  plot <- ggplot(pixelAssocsDT) +
    geom_tile(aes_string(x = "x", y = "y", fill = stat)) +
    geom_path(aes(x=col,y=row,color=area),data = areas, size = 0.5) +
    scale_fill_gradient2() +
    scale_y_reverse() +
    theme_bw() +
    theme(legend.position = "bottom")+
    ggtitle(paste0(trait))

  png(paste0(outDir,"assocsPixelwise_",stat,"_",trait,".png"))
  print(plot)
  dev.off()

}

minpval <- pixelAssocsDT[pVal>0, pVal] %>% min
pixelAssocsDTplot <- pixelAssocsDT %>%
  .[pVal==0, pVal := minpval] %>%
  .[, log10p := (-1)*log10(pVal)]

plot <- ggplot(pixelAssocsDTplot) +
  geom_tile(aes(x = x, y = y, fill = log10p)) +
  geom_path(aes(x=col,y=row,color=area),data = areas,size = 0.5) +
  scale_fill_gradient2() +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "bottom")+
  ggtitle(paste0(trait))

  png(paste0(outDir,"assocsPixelwise_pVal_",trait,".png"))
  print(plot)
  dev.off()

})
