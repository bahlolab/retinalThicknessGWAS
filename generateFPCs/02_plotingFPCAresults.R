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






## Final plots for paper
load("/stornext/Bioinf/data/lab_bahlo/projects/misc/retinalThickness/phenoExploratory/working/cleanedScansFPCA_100fPCs_20221121.RData")
scans <- fread("/stornext/Bioinf/data/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/processedData/scansUnadjustedFinal.csv")
pixels <-  names(scans)[!names(scans) %in% c("patID", "eye", "visit", "sex", "age", "device", "meanRefErr")]
load("/stornext/Bioinf/data/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/rawData/example_data.RData")

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
    ggtitle(paste0("FPC",i)) +
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

##corr plots only, as pdf

lapply(c(1:6), function(i) {

pdf(paste0("/vast/scratch/users/jackson.v/retThickness/fPCs/fpcScoreCorrelations_",i,".pdf"), width = 5.5, height = 6)
(corrPlots[[i]]) %>%
  print
dev.off()

})


## FPC distributions
scoreDT <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/processedData/fpcPhenotypes.txt")

fpcDistPlots <- lapply(c(1:6), function(i) {


  ggplot(scoreDT, aes_string( x = paste0("fpc",i) )) +
    geom_histogram() +
    theme_bw() +
    theme(legend.position = "none")

})

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcScoreDistributions.png", width = 1500, height = 1200)
reduce(fpcDistPlots , `+`) %>%
  print
dev.off()

