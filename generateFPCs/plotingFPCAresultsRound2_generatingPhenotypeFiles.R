library(data.table)
library(magrittr)
library(tidyverse)
library(funData)
library(patchwork)
library(purrr)
library(RColorBrewer)
library(GGally)

load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/cleanedScansFPCAround2_100fPCs_20221124.RData")
scans <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/scansUnadjustedFinalfPCexclusions.csv")
pixels <-  names(scans)[!names(scans) %in% c("patID", "eye", "visit", "sex", "age", "device", "meanRefErr")]

eigen <- data.table(fpc = c(1:100), val=pca$values)
ggplot(eigen, aes(x = fpc, y=val)) +
  geom_line() +
  geom_point()



scoresDT <- pca$scores  %>%
  as.data.table(keep.rownames = T) %>%
  setnames(., c("patID", paste0("fpc",1:100)))

fpcDistPlots <- lapply(c(1:20), function(i) {
  
  
  ggplot(scoresDT, aes_string( x = paste0("fpc",i) )) +
    geom_histogram() +
    theme_bw() +
    theme(legend.position = "none")
  
})

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcRound2ScoreDistributions.png", width = 1500, height = 1200)
reduce(fpcDistPlots , `+`) %>%
  print
dev.off()



load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/scansMatricesUnadjustedFinalfPCexclusions.Data")

## plot scans of outliers


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
    
    png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/outliersRound2_",i,".png"), width = 3000, height = 3000)
    reduce(excludePlots[idx] , `+`) %>%
      print
    dev.off()
  }
  
})


## No exclions needed!!

pcaScoresfilt <- scoresDT

lapply(c(1:25), function(i) {
  
  pc <- paste0("fpc",i)
  top50 <- pcaScoresfilt[order(pcaScoresfilt[, ..pc] , decreasing=TRUE)[1:50], patID]
  bot50 <-pcaScoresfilt[ order(pcaScoresfilt[, ..pc] , decreasing=FALSE)[1:50], patID]
  
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
  
  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpc",i,"Round2_extremes50_outliersRemoved.png"), width = 1000, height = 500)
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
  
  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpc",i,"Round2_top25_outliersRemoved.png"), width = 1250, height = 1250)
  reduce(topPlots25 , `+`) %>%
    print
  dev.off()
  
  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpc",i,"Round2_bot25_outliersRemoved.png"), width = 1250, height = 1250)
  reduce(botPlots25 , `+`) %>%
    print
  dev.off()
  
  
})

fwrite(scoresDT, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/processedData/fpcPhenotypes.txt")