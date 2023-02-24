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

# chr <- 9

sentinels <- lapply(c(1:22), function(chr) {
    
chrSent <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/chr",chr,"sentinelsIDonly_clumpThresh0.001.txt") %>%
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

fwrite(results, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/allSentinelsAllPixelsResults.csv", sep = ",")


  mat <- dcast(results, pixel ~ ID, value.var = "BETA") %>%
    as.matrix(., rownames = "pixel") %>%
    #abs() %>%
    na.omit 
  
  corrMat <- cor(mat)
  
  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/plots/allChr_sentinelsBetasCorrPlot_clumpThresh0.001.png"), width = 2400, height = 2400)
  corrMat %>% corrplot()
  dev.off()
  
  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/plots/allChr_sentinelsPerChr_clumpThresh0.001.png"), width = 800, height = 800)
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
  
  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/plots/allChr_pixelGsentinels.png"), width = 1200, height = 1200)
  print(plot)
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
  
  
  
  # ## N pixels for each SNP
  dt <- results[ , .N, by = POS] %>%
    .[,POS := as.numeric(POS)]
  
  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/plots/fullGrid_chr",chr,"_ChromosomeSigSNPs.png"), width = 1200, height = 600)
  ggplot(dt, aes(x=POS, y = N)) +
    geom_line()
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
    
  
    
    clustersOrdered <- clusterAssignments[,cluster] %>% as.factor
    clustSort <- clusterAssignments[,ID]
    corrMatOrdered <- corrMat[clustSort , clustSort ]
    
    clusterColours <- rep(brewer.pal(8, "Dark2"), 4)
    
    heatmap.2(corrMatOrdered, breaks = seq(0, 1, by = 0.2),
              Rowv = NULL, Colv = NULL, trace = "none",
              ColSideColors = clusterColours[clustersOrdered],
              RowSideColors = clusterColours[clustersOrdered])

    
    
    
    
    
    
    
    
    snpsList <- lapply(sentinels[,ID], function(id) {
      
      
      snpMat <- results[ID == id] %>%
        .[, c("x", "y") := tstrsplit(pixel, "_", type.convert=TRUE)] %>%
        dcast(., y~x, value.var="BETA") %>%
        as.matrix(., rownames="y") %>%
        scale()
      
      return(snpMat)
    })
    names(snpsList) <- sentinels[,ID]
    
    scansFunDataCentred <- snps  - meanFunction(scansFunData)
    
    scansMultiFunData <- multiFunData(list(scansFunDataCentred))
    
    
    snpsArray <- simplify2array(snpsList) %>%
      aperm(., c(3, 2, 1))
    
    y <- dim(snpsArray)[2]
    z <- dim(snpsArray)[3]
    
    domain <- list(c(1:y), c(1:z))
    
    snpsFunData <- funData(domain, snpsArray )
    
    snpsMultiFunData <- multiFunData(list(snpsFunData))
    
   # rm(snpsList, snpsArray, snpsFunData)
    
    pca <- MFPCA(snpsMultiFunData,
                 M = 100,
                 uniExpansions = list(list(type = "splines2Dpen",  k = c(12,12), parallel=F)),
                 verbose = T)
    
    
    
    
    
    
    
    
    snpsListScaled <- lapply(snpsListUnscaled, function(x) {
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
    
    # rm(snpsList, snpsArray, snpsFunData)
    
    pca <- MFPCA(snpsMultiFunData,
                 M = 100,
                 uniExpansions = list(list(type = "splines2Dpen",  k = c(12,12), parallel=F)),
                 verbose = T)
    save(pca, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/allSentinelsFPCA.RData")
    
    
    gap_stat <- clusGap(pca$scores, FUN = kmeans, nstart = 50, K.max = 10, B = 50, spaceH0 = "original")
    k <- maxSE(f = gap_stat$Tab[, "gap"], SE.f = gap_stat$Tab[, "SE.sim"], method = "firstSEmax")
    
    cluster <- kmeans(pca$scores, centers = k)
    
    clusterDT <- data.table(ID = names(cluster$cluster),
                            cluster = cluster$cluster)
    
    
    