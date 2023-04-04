library(data.table)
library(magrittr)
library(funData)
library(MFPCA)
library(dplyr)
library(rlist)
library(cluster)
library(psych)
library(stringr)
library(factoextra)
library(gplots)
library(RColorBrewer)
library(corrplot)
library(ggplot2)
library(doParallel)
library(foreach)
library(patchwork)
library(purrr)


# chr <- 9

sentinels <- lapply(c(1:22), function(chr) {
  
  chrSent <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/chr",chr,"sentinelsIDonly_clumpThresh0.001_withOverlap.txt") %>%
    fread(., header=F)
  return(chrSent)
}) %>%
  rbindlist %>%
  setnames(., "ID")

## read in results for all sentinels
results <- lapply(c(1:22), function(chr) {
  
  dir <- paste0("/vast/scratch/users/jackson.v/retThickness/GWAS/GWsigResults/chr",chr)
  
  print(paste("chromosome",chr))
  chrResults <- lapply(c(1:119), function(slice) {
    
    print(paste(slice))
    
    file <- paste0(dir,"/chr",chr,"Slice",slice,"_5e-5Sig.txt")
    
    sliceResults <- fread(file) %>%
      setnames(., "#POS", "POS") %>%
      .[ID %in% sentinels[,ID]]
    return(sliceResults)
    
  }) %>%
    rbindlist %>%
    .[, chrom := chr]
  
  return(chrResults)
}) %>%
  rbindlist   

fwrite(results, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/allSentinelsAllPixelsResults_clumpThresh0.001_withOverlap.csv", sep = ",")

results <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/allSentinelsAllPixelsResults_clumpThresh0.001_withOverlap.csv")

mat <- dcast(results, pixel ~ ID, value.var = "BETA") %>%
  as.matrix(., rownames = "pixel") %>%
  #abs() %>%
  na.omit 

corrMat <- cor(mat)

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/plots/allChr_sentinelsBetasCorrPlot_clumpThresh0.001_withOverlap.png"), width = 2400, height = 2400)
corrMat %>% corrplot()
dev.off()

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/plots/allChr_sentinelsPerChr_clumpThresh0.001_withOverlap.png"), width = 800, height = 800)
results[, .(ID, chrom)] %>%
  unique %>%
  ggplot(., aes(x=chrom)) +
  geom_bar()
dev.off()

dt <- results[P<5e-8] %>%
  .[ , .N, by = pixel] %>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert = T)] %>%
  .[!is.na(y)]

plot <- ggplot(dt) +
  geom_tile(aes(x = x, y = y, fill = N)) +
  scale_fill_gradient2() +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "bottom")

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/plots/allChr_pixelsGenomeWiseSigSsentinels.png"), width = 1200, height = 1200)
print(plot)
dev.off()

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/plots/allChr_pixels5E-5SigSsentinels.png"), width = 1200, height = 1200)
results[P<5e-5] %>%
  .[ , .N, by = pixel] %>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert = T)] %>%
  .[!is.na(y)] %>%
  ggplot(.) +
  geom_tile(aes(x = x, y = y, fill = N)) +
  scale_fill_gradient2() +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "bottom")
dev.off()

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/plots/allChr_sentinelsPerPixel_clumpThresh0.001.png"), width = 800, height = 800)
ggplot(dt, aes(x=N)) +
  geom_bar()
dev.off()

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/plots/allChr_pixelsPerSentinel_clumpThresh0.001_P5E-8.png"), width = 800, height = 800)
results[P<5e-8] %>%
  .[ , .N, by = ID] %>%
  ggplot(., aes(x=N)) +
  geom_histogram()
dev.off()

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/plots/allChr_pixelsPerSentinel_clumpThresh0.001_P5E-5.png"), width = 800, height = 800)
results[P<5e-5] %>%
  .[ , .N, by = ID] %>%
  ggplot(., aes(x=N)) +
  geom_histogram()
dev.off()







## PCA and cluster on absolute loadings
componentTest <- fa.parallel(corrMat, fa="pc", plot = T, n.obs = nrow(mat))
components <- componentTest$ncomp


pcs <- principal(corrMat, components)

clusters <- factor2cluster(pcs)

clusterDT <- clusters %>%
  as.data.table(., keep.rownames = T) %>%
  melt(., id.vars = "rn")

clusterAssignments <- clusterDT[value!=0]  %>%
  .[, signal := str_remove(variable, "Factor")] %>%
  .[, .(rn, signal)] %>%
  setnames(., c("ID", "cluster"))

pixelsAnnotated <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/annotatedSentinelsPixelwise.txt")

clustersAnnotated <- clusterAssignments[pixelsAnnotated[, 1:18], on = "ID"]



fpcResults <- lapply(c(1:22), function(chr) {
  
  print(paste("chromosome",chr))
  lapply(c(1:8), function(fpc) {
  
  # print(fpc)
  # slice <- 64
  file <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/results/chr",chr,"/chr",chr,"EUR.fpc",fpc,".glm.linear")
  
  fpcResults <- fread(file) %>%
    setnames(., "#CHROM", "CHR") %>%
    .[ID %in% sentinels[,ID]] %>%
    .[, FPC := fpc]
  
  return(fpcResults)
  
}) %>%
  rbindlist
}) %>%
  rbindlist

fwrite(fpcResults, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/FPCresults_allPixelsSentinels.csv", sep = ",")

fpcResults <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/FPCresults_allPixelsSentinels.csv")

fpcLogP <- fpcResults %>%
  .[FPC < 9] %>%
  .[, log10P := (-1)*log(P, 10)] %>%
  .[, sig := case_when(P < 5E-5 ~ "<5E-5",
                       P < 5E-8 ~ "<5E-8",
                       T ~ "")]

# clus <- "RC1"

lapply(paste0("RC", 1:10), function(clus) {
  clusSNPs <- clustersAnnotated[cluster==clus, ID]
  
  snpOrder <- fpcLogP[ID %in% clusSNPs] %>%
    setorder(., CHR, POS) %>%
    .[,ID] %>%
    unique
  
  plot <-  ggplot(fpcLogP[ID %in% clusSNPs], aes(y = ID, x = FPC)) +
    geom_tile(aes(fill = log10P)) +
    scale_fill_gradient2() +
    # scale_y_reverse()  + 
    scale_y_discrete(limits = snpOrder) +
    scale_x_continuous(breaks=seq(1,8,1)) +
    theme_bw() +
    theme(legend.position = "bottom")+
    ggtitle(paste("Cluster",clus)) 
  
  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/plots/FPCresults_cluster",clus,".png"), width = 600, height = 1800)
 plot %>%
    print
  dev.off()
  
})
# geom_text(aes(y=ID,  x = FPC, label=sig)) +
  


clustersOrdered <- clusterAssignments[,cluster] %>% as.factor
clustSort <- clusterAssignments[,ID]
corrMatOrdered <- corrMat[clustSort , clustSort ]

clusterColours <- rep(brewer.pal(8, "Dark2"), 4)

heatmap.2(corrMatOrdered, breaks = seq(0, 1, by = 0.2),
          Rowv = NULL, Colv = NULL, trace = "none",
          ColSideColors = clusterColours[clustersOrdered],
          RowSideColors = clusterColours[clustersOrdered])







nCores <- 6
cluster <- makeCluster(nCores)
doParallel::registerDoParallel(cluster)

# detectedCores <- detectCores()
#
# print(paste("Sense check - there are",detectedCores," cores!"))
set.seed(3467)

## this scale does so by column - not suitable!!

snpsList <- lapply(sentinels[,ID], function(id) {
  
  
  snpMat <- results[ID == id] %>%
    .[, c("x", "y") := tstrsplit(pixel, "_", type.convert=TRUE)] %>%
    dcast(., y~x, value.var="BETA") %>%
    as.matrix(., rownames="y") %>%
    scale()
  
  return(snpMat)
})
names(snpsList) <- sentinels[,ID]

snpsArray <- simplify2array(snpsList) %>%
  aperm(., c(3, 2, 1))

y <- dim(snpsArray)[2]
z <- dim(snpsArray)[3]

domain <- list(c(1:y), c(1:z))

snpsFunData <- funData(domain, snpsArray )

snpsMultiFunData <- multiFunData(list(snpsFunData))
# snpsFunDataCentred <- snpsFunData - meanFunction(snpsFunData)
# 
# snpsMultiFunData <- multiFunData(list(snpsFunDataCentred))

# rm(snpsList, snpsArray, snpsFunData)

pca <- MFPCA(snpsMultiFunData,
             M = 50,
             uniExpansions = list(list(type = "splines2Dpen",  k = c(10,10), parallel=T)),
             verbose = T)










gap_stat <- clusGap(pca$scores, FUN = kmeans, nstart = 50, K.max = 10, B = 50, spaceH0 = "original")
k <- maxSE(f = gap_stat$Tab[, "gap"], SE.f = gap_stat$Tab[, "SE.sim"], method = "firstSEmax")

cluster <- kmeans(pca$scores, centers = k)

clusterDT <- data.table(ID = names(cluster$cluster),
                        cluster = cluster$cluster)






snpsList <- lapply(sentinels[,ID], function(id) {
  
  
  snpMat <- results[ID == id] %>%
    .[, c("x", "y") := tstrsplit(pixel, "_", type.convert=TRUE)] %>%
    dcast(., y~x, value.var="BETA") %>%
    as.matrix(., rownames="y")
  
  return(snpMat)
})
names(snpsList) <- sentinels[,ID]

snpsListScaled <- lapply(snpsList, function(x) {
  m <- mean(x, na.rm=T)
  s <- sd(x, na.rm=T)
  norm <- (x - m) / s
  return(norm)
})

snpsArray <- simplify2array(snpsListScaled) %>%
  aperm(., c(3, 2, 1))

y <- dim(snpsArray)[2]
z <- dim(snpsArray)[3]

domain <- list(c(1:y), c(1:z))

snpsFunData <- funData(domain, snpsArray )

snpsMultiFunData <- multiFunData(list(snpsFunData))


snpsFunDataCentred <- snpsFunData - meanFunction(snpsFunData)

snpsMultiFunData <- multiFunData(list(snpsFunDataCentred))

# rm(snpsList, snpsArray, snpsFunData)

pca <- MFPCA(snpsMultiFunData,
             M = 50,
             uniExpansions = list(list(type = "splines2Dpen",  k = c(10,10), parallel=T)),
             verbose = T)

save(pca, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/allSentinelsFPCA.RData")



mean <- pca$meanFunction[[1]] %>%
  funData::as.data.frame(.) %>%
  as.data.table

ggplot(mean) +
  geom_tile(aes(x = argvals1, y = argvals2, fill = X)) +
  scale_fill_gradient2() +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "bottom")


mask <- mean[,X] %>% is.na %>% which

fpcPlots <- lapply(c(1:8), function(i) {
  
  fpc <- pca$functions[[1]] %>%
    funData::as.data.frame(.) %>%
    as.data.table %>%
    .[obs==i] %>%
    .[!mask] 
  
  plot <- ggplot(fpc) +
    geom_tile(aes(x = argvals1, y = argvals2, fill = X)) +
    scale_fill_gradient2() +
    scale_y_reverse() +
    theme_bw() +
    theme(legend.position = "none")
  # theme(legend.position = "bottom")
  
  return(plot)
  
})

reduce(fpcPlots , `+`) %>%
  print






componentTest <- fa.parallel(corrMat, fa="pc", plot = T, n.obs = nrow(mat))
components <- componentTest$ncomp

pcs <- principal(corrMat, components)

















pcs <- princomp(corrMat, cor=TRUE)
fviz_eig(pcs)

gap_stat <- clusGap(pcs$scores[,1:5], FUN = kmeans, nstart = 50, K.max = 10, B = 50, spaceH0 = "original")
k <- maxSE(f = gap_stat$Tab[, "gap"], SE.f = gap_stat$Tab[, "SE.sim"], method = "firstSEmax")

cluster <- kmeans(pcs$scores, centers = k)

clusterAssignments <- data.table(ID = names(cluster$cluster),
                        cluster = cluster$cluster)

pixelsAnnotated <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/annotatedSentinelsPixelwise.txt")

clustersAnnotated <- clusterAssignments[pixelsAnnotated[, 1:18], on = "ID"]



# fpcResults <- lapply(c(1:22), function(chr) {
#   
#   print(paste("chromosome",chr))
#   lapply(c(1:8), function(fpc) {
#     
#     # print(fpc)
#     # slice <- 64
#     file <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/results/chr",chr,"/chr",chr,"EUR.fpc",fpc,".glm.linear")
#     
#     fpcResults <- fread(file) %>%
#       setnames(., "#CHROM", "CHR") %>%
#       .[ID %in% sentinels[,ID]] %>%
#       .[, FPC := fpc]
#     
#     return(fpcResults)
#     
#   }) %>%
#     rbindlist
# }) %>%
#   rbindlist
# 
# fwrite(fpcResults, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/FPCresults_allPixelsSentinels.csv", sep = ",")

fpcResults <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/FPCresults_allPixelsSentinels.csv")

fpcLogP <- fpcResults %>%
  .[FPC < 9] %>%
  .[, log10P := (-1)*log(P, 10)] %>%
  .[, sig := case_when(P < 5E-5 ~ "<5E-5",
                       P < 5E-8 ~ "<5E-8",
                       T ~ "")] %>%
  .[,logBETA := case_when(BETA > 0 ~ log(BETA+1),
                          BETA < 0 ~ log(abs(BETA)+1)*(-1),
                          T ~ 0)]

# clus <- "RC1"

lapply(c(1:6), function(clus) {
  clusSNPs <- clustersAnnotated[cluster==clus, ID]
  
  snpOrder <- fpcLogP[ID %in% clusSNPs] %>%
    setorder(., CHR, POS) %>%
    .[,ID] %>%
    unique
  
  plot <-  ggplot(fpcLogP[ID %in% clusSNPs], aes(y = ID, x = FPC)) +
    geom_tile(aes(fill = log10P)) +
    scale_fill_gradient2() +
    # scale_y_reverse()  + 
    scale_y_discrete(limits = snpOrder) +
    scale_x_continuous(breaks=seq(1,8,1)) +
    theme_bw() +
    theme(legend.position = "bottom")+
    ggtitle(paste("Cluster",clus)) 
  
  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/plots/FPCresults_cluster",clus,".png"), width = 600, height = 1800)
  plot %>%
    print
  dev.off()
 
  plot <-  ggplot(fpcLogP[ID %in% clusSNPs], aes(y = ID, x = FPC)) +
    geom_tile(aes(fill = logBETA)) +
    scale_fill_gradient2() +
    # scale_y_reverse()  + 
    scale_y_discrete(limits = snpOrder) +
    scale_x_continuous(breaks=seq(1,8,1)) +
    theme_bw() +
    theme(legend.position = "bottom")+
    ggtitle(paste("Cluster",clus)) 
  
  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/plots/FPCresultsBETAs_cluster",clus,".png"), width = 600, height = 1800)
  plot %>%
    print
  dev.off()
  
  
})

fpcLogPclust <- fpcLogP[clustersAnnotated, on = "ID"]

ggplot(fpcLogPclust, aes(BETA, i.BETA, col=as.factor(cluster))) +
  geom_point() +
  facet_wrap(~FPC, scales="free_x")
