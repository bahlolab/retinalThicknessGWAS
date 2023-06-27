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

eigen <- data.table(fpc = c(1:100), 
                    val=pca$values,
                    propExplained=pca$values/sum(pca$values),
                    cumPropExplained=cumsum(pca$values/sum(pca$values))) 
  

screePlots <- lapply(c(100, 50, 30, 15), function (th) {
  
  ggplot(eigen[fpc <= th], aes(x = fpc, y=val)) +
  geom_line() +
  geom_point()
})

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcScree.png", width = 1200, height = 1200)
reduce(screePlots , `+`) %>%
  print
dev.off()


png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcScree_20fPCs.png", width = 600, height = 600)
ggplot(eigen[fpc <= 20], aes(x = fpc, y=val)) +
  geom_line() +
  geom_point() +
  xlab("FPC") +
  ylab("Eigenvalue") + 
  theme_minimal() +
  theme(text = element_text(size = 16))   
dev.off()

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcScree_20fPCs_log.png", width = 600, height = 600)
ggplot(eigen[fpc <= 20], aes(x = fpc, y=val)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(trans = "log10")+
  xlab("FPC") +
  ylab("Eigenvalue") + 
  theme_minimal() +
  theme(text = element_text(size = 16))   


dev.off()

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcCumPropVar_20fPCs.png", width = 600, height = 600)
ggplot(eigen[fpc <= 20], aes(x = fpc, y=cumPropExplained)) +
  geom_line() +
  geom_point()+
  xlab("FPC") +
  ylab("Cumumlative Proportion of Variance Explained") + 
  theme_minimal() +
  theme(text = element_text(size = 16))   

dev.off()

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcCumPropVar_20fPCs_log.png", width = 600, height = 600)
ggplot(eigen[fpc <= 20], aes(x = fpc, y=cumPropExplained)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(trans = "log10")+
  xlab("FPC") +
  ylab("Cumumlative Proportion of Variance Explained") + 
  theme_minimal() +
  theme(text = element_text(size = 16))   

dev.off()

mean <- pca$meanFunction[[1]] %>%
  funData::as.data.frame(.) %>%
  as.data.table %>%
  setnames(., c("obs", "Y", "X", "RT"))
  
png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcMeanFunction.png", width = 600, height = 600)
ggplot(mean) +
  geom_tile(aes(y = Y, x = X, fill = RT)) +
  scale_fill_gradient2() +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "bottom") + 
  labs(fill="Thickness") 
dev.off()



fpcPlots <- lapply(c(1:20), function(i) {

   fpc <- pca$functions[[1]] %>%
    funData::as.data.frame(.) %>%
    as.data.table %>%
    .[obs==i] %>%
     setnames(., c("obs", "Y", "X", "val")) %>%
     .[, pixel := paste(X,Y, sep = "_")] %>%
     .[pixel %in% names(scans)]

  plot <- ggplot(fpc) +
    geom_tile(aes(y = Y, x = X, fill = val)) +
    scale_fill_gradient2() +
    scale_y_reverse() +
    theme_bw() +
    theme(legend.position = "none")
  # theme(legend.position = "bottom")

  return(plot)

})

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcFunctions1-6.png", width = 1200, height = 1800)
((fpcPlots[[1]] + fpcPlots[[2]]) / (fpcPlots[[3]] + fpcPlots[[4]]) / (fpcPlots[[5]] + fpcPlots[[6]])) %>%
  print
dev.off()

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcFunctions1-8.png", width = 1200, height = 2400)
((fpcPlots[[1]] + fpcPlots[[2]]) / (fpcPlots[[3]] + fpcPlots[[4]]) / (fpcPlots[[5]] + fpcPlots[[6]]) / (fpcPlots[[7]] + fpcPlots[[8]])) %>%
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

fpcCorrPlots <- lapply(c(1:8), function(i) {

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

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcScoreCorrelations1-6.png", width = 1200, height = 1800)
((fpcCorrPlots[[1]] + fpcCorrPlots[[2]]) / (fpcCorrPlots[[3]] + fpcCorrPlots[[4]]) / (fpcCorrPlots[[5]] + fpcCorrPlots[[6]])) %>%
  print
dev.off()

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcScoreCorrelations1-8.png", width = 1200, height = 2400)
((fpcCorrPlots[[1]] + fpcCorrPlots[[2]]) / (fpcCorrPlots[[3]] + fpcCorrPlots[[4]]) / (fpcCorrPlots[[5]] + fpcCorrPlots[[6]]) / (fpcCorrPlots[[7]] + fpcCorrPlots[[8]])) %>%
  print
dev.off()





scoreDT <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/processedData/fpcPhenotypes.txt")

fpcDistPlots <- lapply(c(1:8), function(i) {


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



scoresDT <- pca$scores  %>%
  as.data.table(keep.rownames = T) %>%
  setnames(., c("patID", paste0("fpc",1:100)))

## plot scans of outliers

# exclude <- lapply(c(1:100), function(i){
#   
#   col <- paste0("fpc",i)
#   vals <- scoresDT[, get(col)] %>% as.vector
#   
#   m <- vals  %>% mean
#   sd <- vals  %>% sd
#   
#   if(i <= 6) {
#     exclude <- ifelse(vals %between% c((m-(8*sd)), (m+(8*sd))), 0, 1)
#   } else{
#     exclude <- ifelse(vals %between% c((m-(5*sd)), (m+(5*sd))), 0, 1)
#   }
#   
# })
# 
# excludeAll <- list.cbind(exclude) %>%
#   as.data.table %>%
#   cbind(scoresDT[, "patID"], .) %>%
#   .[, excSum := rowSums(.SD, na.rm=T), .SDcols = paste0("V",c(1:100))] %>%
#   .[, exclude := ifelse(excSum > 0, 1, 0)]
# 
# exclIDs <- excludeAll[exclude==1, patID]
# 
# 
# 
# excludeList <- scansList[exclIDs]
# 
# excludePlots <- lapply(excludeList, function(scan) {
#   
#   scanDT <- scan %>%
#     melt %>%
#     as.data.table
#   scanPlot <-  ggplot(scanDT) +
#     geom_tile(aes(x = Var2, y = Var1, fill = value)) +
#     scale_fill_continuous(limits = c(40, 140), 
#                           low = brewer.pal(6, "Blues")[1], 
#                           high = brewer.pal(6, "Blues")[6]) +
#     scale_y_reverse() +
#     theme_bw() +
#     theme(legend.position = "bottom")
#   
#   return(scanPlot)
# })

# lapply(c(1:100), function(i) {
#   
# pc <- paste0("V",i)
# idx <- excludeAll[get(pc)==1 & exclude==1, patID]
# 
# if(length(idx) > 0){
#   
# png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/outliers_",i,".png"), width = 3000, height = 3000)
# reduce(excludePlots[idx] , `+`) %>%
#   print
# dev.off()
# }
# 
# })
# 

# pcaScoresfilt <- pca$scores %>%
#   as.data.table(., keep.rownames = T) %>%
#   .[which(!rn %in% exclIDs)]

pcaScoresfilt <-scoresDT

lapply(c(1:8), function(i) {
  
  pc <- paste0("fpc",i)
  top100 <- pcaScoresfilt[order(pcaScoresfilt[, ..pc] , decreasing=TRUE)[1:100], patID]
  bot100 <-pcaScoresfilt[ order(pcaScoresfilt[, ..pc] , decreasing=FALSE)[1:100], patID]
  
  top100List <- scansList[top100]
  top100Mean <- Reduce("+", top100List) / length(top100List) 
  
  bot100List <- scansList[bot100]
  bot100Mean <- Reduce("+", bot100List) / length(bot100List) 
  
  topDT <- top100Mean %>%
    melt %>%
    as.data.table
  
  botDT <- bot100Mean %>%
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
  
  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpc",i,"_extremes100.png"), width = 1000, height = 500)
  print(topPlot + botPlot)
  dev.off()
  
  
  # topPlots25 <- lapply(top100List[1:25], function(scan) {
  #   
  #   scanDT <- scan %>%
  #     melt %>%
  #     as.data.table
  #   scanPlot <-  ggplot(scanDT) +
  #     geom_tile(aes(x = Var2, y = Var1, fill = value)) +
  #     scale_fill_continuous(limits = c(40, 140), 
  #                           low = brewer.pal(6, "Blues")[1], 
  #                           high = brewer.pal(6, "Blues")[6]) +
  #     scale_y_reverse() +
  #     theme_bw() +
  #     theme(legend.position = "bottom")
  #   
  #   return(scanPlot)
  # })
  # 
  # botPlots25 <- lapply(bot100List[1:25], function(scan) {
  #   
  #   scanDT <- scan %>%
  #     melt %>%
  #     as.data.table
  #   scanPlot <-  ggplot(scanDT) +
  #     geom_tile(aes(x = Var2, y = Var1, fill = value)) +
  #     scale_fill_continuous(limits = c(40, 140), 
  #                           low = brewer.pal(6, "Blues")[1], 
  #                           high = brewer.pal(6, "Blues")[6]) +
  #     scale_y_reverse() +
  #     theme_bw() +
  #     theme(legend.position = "bottom")
  #   
  #   return(scanPlot)
  # })
  # 
  # png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpc",i,"_top25_outliersRemoved.png"), width = 1250, height = 1250)
  # reduce(topPlots25 , `+`) %>%
  #   print
  # dev.off()
  # 
  # png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpc",i,"_bot25_outliersRemoved.png"), width = 1250, height = 1250)
  # reduce(botPlots25 , `+`) %>%
  #   print
  # dev.off()
  # 
  
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



## Final plots for paper
load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenoExploratory/working/cleanedScansFPCA_100fPCs_20221121.RData")
scans <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/processedData/scansUnadjustedFinal.csv")
pixels <-  names(scans)[!names(scans) %in% c("patID", "eye", "visit", "sex", "age", "device", "meanRefErr")]
load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/rawData/example_data.RData")

pcaScoresfilt <- pca$scores  %>%
  as.data.table(keep.rownames = T) %>%
  setnames(., c("patID", paste0("fpc",1:100)))
rm(pca)

corrPlots <- lapply(c(1:6), function(i) {

  pc <- paste0("fpc",i)
  print(paste(pc))

fpcCorr <- cor(scans[, ..pixels], pcaScoresfilt[, ..pc]) %>%
  as.data.table(keep.rownames = T) %>%
  .[, c("y", "x") := tstrsplit(rn, "_", type.convert=TRUE)]

setnames(fpcCorr, pc, "cor")

  corrPlot <- ggplot(fpcCorr) +
  geom_tile(aes(x = x, y = y, fill = cor)) +
    scale_fill_gradient2() +
    scale_y_reverse() +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 16))+
    ggtitle("(i)") +
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

  return(corrPlot)
})


rm(scans)
load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/processedData/scansMatricesUnadjustedFinal.RData")

extremePlots <- lapply(c(1:6), function(i) {

  pc <- paste0("fpc",i)
  print(paste(pc))

  top100 <- pcaScoresfilt[order(pcaScoresfilt[, ..pc] , decreasing=TRUE)[1:100], patID]
  bot100 <-pcaScoresfilt[ order(pcaScoresfilt[, ..pc] , decreasing=FALSE)[1:100], patID]
  
  top100List <- scansList[top100]
  top100Mean <- Reduce("+", top100List) / length(top100List) 
  
  bot100List <- scansList[bot100]
  bot100Mean <- Reduce("+", bot100List) / length(bot100List) 
  
  topDT <- top100Mean %>%
    melt %>%
    as.data.table %>%
    setnames(., c("y", "x", "RT")) %>%
    .[, subSet := "(ii)"]
  
  botDT <- bot100Mean %>%
    melt %>%
    as.data.table%>%
    setnames(., c("y", "x", "RT")) %>%
    .[, subSet := "(iii)"]
  
  extremes <- rbind(topDT, botDT)

  minThick <- min(extremes[,RT], na.rm=T) %>% floor(.)
  maxThick <- max(extremes[,RT], na.rm=T) %>% ceiling(.)
  
  extremesPlot <-  ggplot(extremes) +
    geom_tile(aes(x = x, y = y, fill = RT)) +
     facet_wrap(vars(subSet), nrow = 2) +                      
     scale_fill_continuous(limits = c(minThick, maxThick), 
                          low = brewer.pal(9, "Blues")[1], 
                          high = brewer.pal(9, "Blues")[9]) +                           
    scale_y_reverse() +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 16),
    strip.background = element_blank(), 
    strip.text.x = element_text(size = 16, hjust = 0)) +
    ggtitle("") +
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

return(extremesPlot)
})

lapply(c(1:6), function(i) {

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcScoreCorrelationsTopBot_",i,".png"), width = 900, height = 650)
(corrPlots[[i]] + extremePlots[[i]] + 
  plot_layout(widths = c(2, 1))) %>%
  print
dev.off()

})
