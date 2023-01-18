library(data.table)
library(magrittr)
library(tidyverse)
library(funData)
library(patchwork)
library(purrr)
library(RColorBrewer)
library(GGally)

# load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenoExploratory/working/cleanFinalScans.RData")
load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenoExploratory/working/cleanedScansFPCA_100fPCs_20221121.RData")
scans <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/processedData/scansUnadjustedFinal.csv")
pixels <-  names(scans)[!names(scans) %in% c("patID", "eye", "visit", "sex", "age", "device", "meanRefErr")]

eigen <- data.table(fpc = c(1:100), val=pca$values)

screePlots <- lapply(c(100, 50, 30, 15), function (th) {
  
  ggplot(eigen[fpc <= th], aes(x = fpc, y=val)) +
  geom_line() +
  geom_point()
})

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcScree.png", width = 1200, height = 1200)
reduce(screePlots , `+`) %>%
  print
dev.off()

mean <- pca$meanFunction[[1]] %>%
  funData::as.data.frame(.) %>%
  as.data.table
  
png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcMeanFunction.png", width = 600, height = 600)
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

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcFunctions.png", width = 1500, height = 1200)
reduce(fpcPlots , `+`) %>%
  print
dev.off()




# ret <- scansList[[4886]] %>%
#   melt %>%
#   as.data.table

# png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/exampleRetina.png", width = 1500, height = 1200)
# ggplot(ret) +
#   geom_tile(aes(x = Var1, y = Var2, fill = value)) +
#   scale_fill_gradient2() +
#   scale_y_reverse() +
#   theme_bw() +
#   theme(legend.position = "bottom")
# dev.off()

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

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcScoreCorrelations.png", width = 1500, height = 1200)
reduce(fpcCorrPlots , `+`) %>%
print
dev.off()


scoreDT <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/processedData/fpcPhenotypes.txt")

fpcDistPlots <- lapply(c(1:20), function(i) {


  ggplot(scoreDT, aes_string( x = paste0("fpc",i) )) +
    geom_histogram() +
    theme_bw() +
    theme(legend.position = "none")

})

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcScoreDistributions.png", width = 1500, height = 1200)
reduce(fpcDistPlots , `+`) %>%
  print
dev.off()



load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/processedData/scansMatricesUnadjustedFinal.RData")



## plot scans of outliers
scoresDT <- pca$scores  %>%
  as.data.table(keep.rownames = T) %>%
  setnames(., c("patID", paste0("fpc",1:100)))


exclude <- lapply(c(1:100), function(i){
  
  col <- paste0("fpc",i)
  vals <- scoresDT[, get(col)] %>% as.vector
  
  m <- vals  %>% mean
  sd <- vals  %>% sd
  
  if(i <= 6) {
    exclude <- ifelse(vals %between% c((m-(8*sd)), (m+(8*sd))), 0, 1)
  } else{
    exclude <- ifelse(vals %between% c((m-(5*sd)), (m+(5*sd))), 0, 1)
  }
  
})

excludeAll <- list.cbind(exclude) %>%
  as.data.table %>%
  cbind(scoresDT[, "patID"], .) %>%
  .[, excSum := rowSums(.SD, na.rm=T), .SDcols = paste0("V",c(1:100))] %>%
  .[, exclude := ifelse(excSum > 0, 1, 0)]

exclIDs <- excludeAll[exclude==1, patID]



excludeList <- scansList[exclIDs]

excludePlots <- lapply(excludeList, function(scan) {
  
  scanDT <- scan %>%
    melt %>%
    as.data.table
  scanPlot <-  ggplot(scanDT) +
    geom_tile(aes(x = Var2, y = Var1, fill = value)) +
    scale_fill_continuous(limits = c(40, 140), 
                          low = brewer.pal(6, "Blues")[1], 
                          high = brewer.pal(6, "Blues")[6]) +
    scale_y_reverse() +
    theme_bw() +
    theme(legend.position = "bottom")
  
  return(scanPlot)
})

lapply(c(1:100), function(i) {
  
pc <- paste0("V",i)
idx <- excludeAll[get(pc)==1 & exclude==1, patID]

if(length(idx) > 0){
  
png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/outliers_",i,".png"), width = 3000, height = 3000)
reduce(excludePlots[idx] , `+`) %>%
  print
dev.off()
}

})


pcaScoresfilt <- pca$scores %>%
  as.data.table(., keep.rownames = T) %>%
  .[which(!rn %in% exclIDs)]

lapply(c(1:25), function(i) {
  
  pc <- paste0("V",i)
  top50 <- pcaScoresfilt[order(pcaScoresfilt[, ..pc] , decreasing=TRUE)[1:50], rn]
  bot50 <-pcaScoresfilt[ order(pcaScoresfilt[, ..pc] , decreasing=FALSE)[1:50], rn]
  
  top50List <- scansList[top50]
  top50Mean <- Reduce("+", top50List) / length(top50List) 
  
  bot50List <- scansList[bot50]
  bot50Mean <- Reduce("+", bot50List) / length(bot50List) 
  
  topDT <- top50Mean %>%
    melt %>%
    as.data.table
  
  botDT <- bot50Mean %>%
    melt %>%
    as.data.table
  
  minThick <- min(min(topDT[,value], na.rm=T), min(botDT[,value], na.rm=T), na.rm=T) %>% floor(.)
  maxThick <- max(max(topDT[,value], na.rm=T), max(botDT[,value], na.rm=T), na.rm=T) %>% ceiling(.)
  
  topPlot <-  ggplot(topDT) +
    geom_tile(aes(x = Var2, y = Var1, fill = value)) +
    scale_fill_continuous(limits = c(minThick, maxThick), 
                          low = brewer.pal(6, "Blues")[1], 
                          high = brewer.pal(6, "Blues")[6]) +
    scale_y_reverse() +
    theme_bw() +
    theme(legend.position = "bottom")
  
  botPlot <-  ggplot(botDT) +
    geom_tile(aes(x = Var2, y = Var1, fill = value)) +
    scale_fill_continuous(limits = c(minThick, maxThick), 
                          low = brewer.pal(6, "Blues")[1], 
                          high = brewer.pal(6, "Blues")[6]) +
    scale_y_reverse() +
    theme_bw() +
    theme(legend.position = "bottom")
  
  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpc",i,"_extremes50_outliersRemoved.png"), width = 1000, height = 500)
  print(topPlot + botPlot)
  dev.off()
  
  
  topPlots25 <- lapply(top50List[1:25], function(scan) {
    
    scanDT <- scan %>%
      melt %>%
      as.data.table
    scanPlot <-  ggplot(scanDT) +
      geom_tile(aes(x = Var2, y = Var1, fill = value)) +
      scale_fill_continuous(limits = c(40, 140), 
                            low = brewer.pal(6, "Blues")[1], 
                            high = brewer.pal(6, "Blues")[6]) +
      scale_y_reverse() +
      theme_bw() +
      theme(legend.position = "bottom")
    
    return(scanPlot)
  })
  
  botPlots25 <- lapply(bot50List[1:25], function(scan) {
    
    scanDT <- scan %>%
      melt %>%
      as.data.table
    scanPlot <-  ggplot(scanDT) +
      geom_tile(aes(x = Var2, y = Var1, fill = value)) +
      scale_fill_continuous(limits = c(40, 140), 
                            low = brewer.pal(6, "Blues")[1], 
                            high = brewer.pal(6, "Blues")[6]) +
      scale_y_reverse() +
      theme_bw() +
      theme(legend.position = "bottom")
    
    return(scanPlot)
  })
  
  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpc",i,"_top25_outliersRemoved.png"), width = 1250, height = 1250)
  reduce(topPlots25 , `+`) %>%
    print
  dev.off()
  
  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpc",i,"_bot25_outliersRemoved.png"), width = 1250, height = 1250)
  reduce(botPlots25 , `+`) %>%
    print
  dev.off()
  
  
})







scoreDT <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fPCscores_noExclusions.csv")

fpcDistPlots <- lapply(c(1:25), function(i) {
  
  
  ggplot(scoreDT, aes_string( x = paste0("fpc",i) )) +
    geom_histogram() +
    theme_bw() +
    theme(legend.position = "none")
  
})

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcScoreDistributionsNoExclusions.png", width = 1500, height = 1500)
reduce(fpcDistPlots , `+`) %>%
  print
dev.off()


scoreDT <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/processedData/fpcPhenotypes.txt")

fpcDistPlots <- lapply(c(1:25), function(i) {
  
  
  ggplot(scoreDT, aes_string( x = paste0("fpc",i) )) +
    geom_histogram() +
    theme_bw() +
    theme(legend.position = "none")
  
})

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcScoreDistributions.png", width = 1500, height = 1500)
reduce(fpcDistPlots , `+`) %>%
  print
dev.off()
