#!/usr/bin/env Rscript

## Run using R/4.1.2

library(data.table)
library(magrittr)
library(dplyr)
library(parallel)
library(rlist)
library(ggplot2)
library(DescTools)


scansDT  <-fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/processedData/scansUnadjustedFinal.csv")

pixels <- names(scansDT)[!names(scansDT) %in% c("patID", "visit", "eye", "sex", "age", "device", "meanRefErr")]


models <- lapply(pixels, function(y) {

  pixelIdx <- which(pixels==y)
  paste("pixel", pixelIdx) %>% print

  model <- lm(get(y) ~ meanRefErr, data = scansDT, na.action=na.exclude)
  # return(model)

  beta <- model$coefficients["meanRefErr"]
  t <- summary(model)$coefficients["meanRefErr", "t value"]
  pVal <- summary(model)$coefficients["meanRefErr", "Pr(>|t|)"]

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

save(models, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/reffErrorModels_alt.RData")


pixelAssocsDT <- list.map(models, modelStats) %>%
  rbindlist %>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert=TRUE)]


for(stat in c("beta", "t")) {

  plot <- ggplot(pixelAssocsDT , aes_string(x = "x", y = "y", fill = stat)) +
    geom_tile() +
    scale_fill_gradient2()

  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/plots/refErrorPixelwise_",stat,"_alt.png"))
  print(plot)
  dev.off()

}

minpval <- pixelAssocsDT[pVal>0, pVal] %>% min
pixelAssocsDTplot <- pixelAssocsDT %>%
  .[pVal==0, pVal := minpval] %>%
  .[, log10p := (-1)*log10(pVal)]


png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/plots/refErrorPixelwise_pVal_alt.png")
ggplot(pixelAssocsDTplot , aes(x = x, y = y, fill = log10p)) +
  geom_tile()
dev.off()
