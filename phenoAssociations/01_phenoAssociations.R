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
# load(paste0(dataDir,"Data_for_retinal_areas.RData"))

pixels <- names(scansDT)[!names(scansDT) %in% c("patID", "visit", "eye", "sex", "age", "device", "meanRefErr")]


 linkage <- fread(paste0(dataDir,"idLinkage.txt"))
 scansDTlinked <- linkage[scansDT, on = c("patID" = "patID")]

file <- paste0(dataDir,"ukb41258.tab")
names <- fread(file, nrows = 0) %>%
  names

getIdx <- function(code) {
    names %like any% code %>% which
  }

cols <- getIdx(c("f.eid", "f.48.0.0", "f.48.1.0", "f.49.0.0", "f.49.1.0",
"f.50.0.0", "f.50.1.0", "f.51.0.0", "f.51.1.0", "f.93.0.%", "f.93.1.%",
"f.94.0.%", "f.94.1.%", "f.4079.0.%", "f.4079.1.%", "f.4080.0.%", "f.4080.1.%",
 "f.6177.0.%", "f.6177.1.%", "f.21002.%", "f.21001.%", "f.20002.%"))

covs <- fread(file, select = cols)

# pcs <- paste0(dataDir,"ukb32825.tab") %>%
#   fread(., select = c(1, 1145:1154)) %>%
#   setnames(. ,c("patIDhda", paste0("PC",c(1:10))))

ancestries <- lapply(c("EUR", "CSA", "AFR"), function(anc) {
  
  ids <- paste0(dataDir,"sampleList_singleIDs_",anc,".txt") %>%
    fread %>%
    .[, ancestry := anc] 
  
return(ids)
  
}) %>%
  rbindlist %>%
  .[, ancestry := as.factor(ancestry)] %>%
  .[, ancestry := relevel(ancestry, ref = "EUR")]

covsMerged <- scansDTlinked[, .(patID, patIDhda, visit)] %>%
  ancestries[., on = c(IID = "patIDhda")] %>%
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
  .[, MS := ifelse(rowSums(.SD == "1261", na.rm = TRUE) > 0, 1, 0), .SDcols = patterns("^f\\.20002\\..*\\..*$")] %>%
  setnames(., "f.eid", "patID")

phenoOut <- covsMerged[, .(patID, ancestry, visit, waist, hip, standHeight, sitHeight, weight, bmi, rawSBP, rawDBP, adjustedSBP, adjustedDBP, BPmedication, MS)] %>%
  .[scansDTlinked, on = c("patID", "visit")] %>%
  .[, ageSq := age^2] %>%
  .[!is.na(ancestry)] 


## tidy up
rm(scansDT)
rm(scansDTlinked)

## summarise demographics
phenos <- c("sex", "meanRefErr", "age", "ageSq", "standHeight", "weight", "bmi", "MS")

    col_means <- phenoOut[, lapply(.SD, mean, na.rm = T), by = ancestry, .SDcols = c("sex", "meanRefErr", "age",  "standHeight", "weight", "bmi")]
    col_sds <- phenoOut[, lapply(.SD, sd, na.rm = T), by = ancestry, .SDcols = c("meanRefErr", "age",  "standHeight", "weight", "bmi")]

result <- col_means[col_sds, on = "ancestry"] %>%
setnames(.,  c("ancestry", "propMale", 
paste0( c("meanRefErr", "age",  "standHeight", "weight", "bmi"), "Mean"), 
paste0( c("meanRefErr", "age",  "standHeight", "weight", "bmi"), "SD")))

# Print the result
print(result)

## Overall mean and var per pixel
pixMeans <- phenoOut[, lapply(.SD, mean, na.rm = T),  .SDcols = pixels]
pixVars <- phenoOut[, lapply(.SD, var, na.rm = T),  .SDcols = pixels]

pixSumms <- data.table(pixel = pixels,
mean = pixMeans %>% unlist,
var = pixVars %>% unlist) %>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert=TRUE)] 

load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/rawData/example_data.RData")

for(stat in c("mean", "var")) {
  
  plot <- ggplot(pixSumms) +
    geom_tile(aes_string(x = "x", y = "y", fill = stat)) +
    scale_fill_gradient2() +
    scale_y_reverse() +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 16))+
    ggtitle(paste0(stat)) +
      geom_path(aes(x=col,y=diag1),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=diag2),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=line1),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=line2),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=row,group=area),color="grey50",data = areas[areas$edtrs,],size = 0.5) +
     guides(fill = guide_colorbar(label.position = "bottom",
                               title.position = "left", 
                               # draw border around the legend
                               frame.colour = "black",
                               barwidth = 15,
                               barheight = 1.5)) 

  png(paste0(outDir,"summaryPixelwise_",stat,"_overall.png"), height = 650, width = 600)
  print(plot)
  dev.off()

}


##  mean and var per pixel, by sex
pixMeans <- phenoOut[, lapply(.SD, mean, na.rm = T), by = sex, .SDcols = pixels] %>%
 melt(., measure.vars = pixels, value.name = "mean", variable.name = "pixel")%>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert=TRUE)] 

pixVars <- phenoOut[, lapply(.SD, var, na.rm = T), by = sex, .SDcols = pixels]  %>%
 melt(., measure.vars = pixels, value.name = "var", variable.name = "pixel") %>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert=TRUE)] 


minMean <- pixMeans[,mean] %>% min(., na.rm=T) %>% floor
minVar <- pixVars[,var] %>% min(., na.rm=T) %>% floor
maxMean <- pixMeans[,mean] %>% max(., na.rm=T) %>% ceiling
maxVar <- pixVars[,var] %>% max(., na.rm=T) %>% ceiling

load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/rawData/example_data.RData")

for(val in c(0,1)) {
  
  text <- ifelse(val==1, "Male", "Female")

  plot <- ggplot(pixMeans[sex==val]) +
    geom_tile(aes_string(x = "x", y = "y", fill = "mean")) +
    scale_fill_gradient2(limits = c(minMean, maxMean)) +
    scale_y_reverse() +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 16))+
    ggtitle(paste0(text, " - Mean")) +
      geom_path(aes(x=col,y=diag1),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=diag2),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=line1),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=line2),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=row,group=area),color="grey50",data = areas[areas$edtrs,],size = 0.5) +
     guides(fill = guide_colorbar(label.position = "bottom",
                               title.position = "left", 
                               # draw border around the legend
                               frame.colour = "black",
                               barwidth = 15,
                               barheight = 1.5)) 

  png(paste0(outDir,"summaryPixelwise_mean_",text,".png"), height = 650, width = 600)
  print(plot)
  dev.off()

  plot <- ggplot(pixVars[sex==val]) +
    geom_tile(aes_string(x = "x", y = "y", fill = "var")) +
    scale_fill_gradient2(limits = c(minVar, maxVar)) +
    scale_y_reverse() +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 16))+
    ggtitle(paste0(text, " - Var")) +
      geom_path(aes(x=col,y=diag1),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=diag2),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=line1),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=line2),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=row,group=area),color="grey50",data = areas[areas$edtrs,],size = 0.5) +
     guides(fill = guide_colorbar(label.position = "bottom",
                               title.position = "left", 
                               # draw border around the legend
                               frame.colour = "black",
                               barwidth = 15,
                               barheight = 1.5)) 

  png(paste0(outDir,"summaryPixelwise_vars_",text,".png"), height = 650, width = 600)
  print(plot)
  dev.off()

}






##  mean and var per pixel, by age category

phenoOut <- phenoOut[ ,ageCat := case_when(age < 45 ~ "<45",
age >= 45 & age < 55 ~ "45-54",
age >= 55 & age < 65 ~ "55-64",
age >= 65 ~ "65+",
T ~ NA_character_)]


pixMeans <- phenoOut[, lapply(.SD, mean, na.rm = T), by = ageCat, .SDcols = pixels] %>%
 melt(., measure.vars = pixels, value.name = "mean", variable.name = "pixel")%>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert=TRUE)] 

pixVars <- phenoOut[, lapply(.SD, var, na.rm = T), by = ageCat, .SDcols = pixels]  %>%
 melt(., measure.vars = pixels, value.name = "var", variable.name = "pixel") %>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert=TRUE)] 


minMean <- pixMeans[,mean] %>% min(., na.rm=T) %>% floor
minVar <- pixVars[,var] %>% min(., na.rm=T) %>% floor
maxMean <- pixMeans[,mean] %>% max(., na.rm=T) %>% ceiling
maxVar <- pixVars[,var] %>% max(., na.rm=T) %>% ceiling

load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/rawData/example_data.RData")

for(val in phenoOut[,ageCat] %>% unique) {
  

  plot <- ggplot(pixMeans[ageCat==val]) +
    geom_tile(aes_string(x = "x", y = "y", fill = "mean")) +
    scale_fill_gradient2(limits = c(minMean, maxMean)) +
    scale_y_reverse() +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 16))+
    ggtitle(paste0(val, " - Mean")) +
      geom_path(aes(x=col,y=diag1),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=diag2),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=line1),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=line2),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=row,group=area),color="grey50",data = areas[areas$edtrs,],size = 0.5) +
     guides(fill = guide_colorbar(label.position = "bottom",
                               title.position = "left", 
                               # draw border around the legend
                               frame.colour = "black",
                               barwidth = 15,
                               barheight = 1.5)) 

  png(paste0(outDir,"summaryPixelwise_mean_",val,".png"), height = 650, width = 600)
  print(plot)
  dev.off()

  plot <- ggplot(pixVars[ageCat==val]) +
    geom_tile(aes_string(x = "x", y = "y", fill = "var")) +
    scale_fill_gradient2(limits = c(minVar, maxVar)) +
    scale_y_reverse() +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 16))+
    ggtitle(paste0(val, " - Var")) +
      geom_path(aes(x=col,y=diag1),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=diag2),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=line1),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=line2),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=row,group=area),color="grey50",data = areas[areas$edtrs,],size = 0.5) +
     guides(fill = guide_colorbar(label.position = "bottom",
                               title.position = "left", 
                               # draw border around the legend
                               frame.colour = "black",
                               barwidth = 15,
                               barheight = 1.5)) 

  png(paste0(outDir,"summaryPixelwise_vars_",val,".png"), height = 650, width = 600)
  print(plot)
  dev.off()

}


## mean by both sex and age
pixMeans <- phenoOut[, lapply(.SD, mean, na.rm = T), by = c("sex", "ageCat"), .SDcols = pixels] %>%
 melt(., measure.vars = pixels, value.name = "mean", variable.name = "pixel")%>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert=TRUE)] 


minMean <- pixMeans[,mean] %>% min(., na.rm=T) %>% floor
maxMean <- pixMeans[,mean] %>% max(., na.rm=T) %>% ceiling

load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/rawData/example_data.RData")

for(val in phenoOut[,ageCat] %>% unique) {
for(sexVal in c(0,1)) {
  
  text <- ifelse(sexVal==1, "Male", "Female")

 

  plot <- ggplot(pixMeans[ageCat==val & sex==sexVal]) +
    geom_tile(aes_string(x = "x", y = "y", fill = "mean")) +
    scale_fill_gradient2(limits = c(minMean, maxMean)) +
    scale_y_reverse() +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 16))+
    ggtitle(paste0(text," - ",val, " - Mean")) +
      geom_path(aes(x=col,y=diag1),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=diag2),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=line1),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=line2),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=row,group=area),color="grey50",data = areas[areas$edtrs,],size = 0.5) +
     guides(fill = guide_colorbar(label.position = "bottom",
                               title.position = "left", 
                               # draw border around the legend
                               frame.colour = "black",
                               barwidth = 15,
                               barheight = 1.5)) 

  png(paste0(outDir,"summaryPixelwise_mean_",val,"_",text,".png"), height = 650, width = 600)
  print(plot)
  dev.off()


}
}






# phenos <- c("sex",  "age", "ageSq", "waist", "hip", "standHeight", "sitHeight", "weight", "bmi", "rawSBP", "rawDBP", "adjustedSBP", "adjustedDBP", "BPmedication", "MS")



vals <- c("ancestryAFR", "ancestryCSA", "sex", "meanRefErr", "age", "ageSq", "standHeight")

assocs <-  lapply(pixels, function(y) {

  pixelIdx <- which(pixels==y)
  paste("pixel", pixelIdx) %>% print

  model <- lm(get(y) ~ ancestry + sex + meanRefErr + age + ageSq + standHeight + device + eye, data = phenoOut, na.action=na.exclude)

  betas <- model$coefficients[vals]
  pVals <- summary(model)$coefficients[vals, "Pr(>|t|)"]

  out<- data.table(pixel = y,
  variable =  vals,
    beta = betas,
    pVal = pVals)

return(out)

}) %>%
rbindlist


minpval <- assocs[pVal>0, pVal] %>% min

assocs  <- assocs %>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert=TRUE)] %>%
  .[pVal==0, pVal := minpval] %>%
  .[, log10p := (-1)*log10(pVal)]

load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/rawData/example_data.RData")

for(stat in c("beta", "log10p")) {
  
  lapply(vals, function(value) {

  plot <- ggplot(assocs[variable==value]) +
    geom_tile(aes_string(x = "x", y = "y", fill = stat)) +
    scale_fill_gradient2() +
    scale_y_reverse() +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 16))+
    ggtitle(paste0(value)) +
      geom_path(aes(x=col,y=diag1),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=diag2),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=line1),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=line2),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=row,group=area),color="grey50",data = areas[areas$edtrs,],size = 0.5) +
     guides(fill = guide_colorbar(label.position = "bottom",
                               title.position = "left", 
                               # draw border around the legend
                               frame.colour = "black",
                               barwidth = 15,
                               barheight = 1.5)) 

  png(paste0(outDir,"assocsPixelwise_",stat,"_",value,".png"), height = 650, width = 600)
  print(plot)
  dev.off()

})
}


























