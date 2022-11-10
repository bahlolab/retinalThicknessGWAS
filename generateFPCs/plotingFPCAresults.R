library(data.table)
library(magrittr)
library(tidyverse)
library(funData)
library(patchwork)
library(purrr)

# load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenoExploratory/working/cleanFinalScans.RData")
load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenoExploratory/working/cleanedScansFPCA_20220927.RData")
scans <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/processedData/scansUnadjustedFinal.csv")
pixels <-  names(scans)[!names(scans) %in% c("patID", "eye", "visit", "sex", "age", "device", "meanRefErr")]

eigen <- data.table(fpc = c(1:100), val=pca$values)
ggplot(eigen, aes(x = fpc, y=val)) +
  geom_line() +
  geom_point()

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
