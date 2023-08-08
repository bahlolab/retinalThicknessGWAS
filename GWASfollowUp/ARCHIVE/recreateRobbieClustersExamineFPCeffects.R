## pixel clusters

load(file="~/lab_bahlo/projects/mactel/PRJ2022002/analysis_versions/version001/Post-GWAS/processed_data/Cleaned_SNPs_results.RData")

sce_snp$Biomarker <- sce_snp$ID

pixel_list <- rowData(sce_snp) %>% as.data.frame()

sce_snp

logcounts <- logcounts(sce_snp)

mat <- matrix(1,ncol=2,nrow=2)
t(t(mat)/c(1,2))

logcounts <- t(t(logcounts)/ colMaxs(abs(logcounts)))
logcounts(sce_snp) <- logcounts
logcounts(sce_pixels) <- t(logcounts)


library(dendextend)
# Check if HClust can do a better job
mat <- logcounts(sce_snp)
mat <- t(mat)
hclust <- hclust(dist(mat),method = "ward.D2")
dend <-as.dendrogram(hclust)
k <- 10
clust.cutree <- dendextend:::cutree(dend, k=k, order_clusters_as_data = FALSE)


# Assign the Hcluster to the dataset
clusts <- clust.cutree
clusts <- data.frame(Biomarker=names(clusts),cluster_hclust=clusts)
clusters_dat <- colData(sce_snp) %>% as.data.frame()
clusters_dat <- left_join(clusters_dat,clusts)
#table(clusters_dat$cluster_snp,clusters_dat$cluster_hclust)
sce_snp$cluster_hclust <- as.character(clusters_dat$cluster_hclust)

## output clusters
pixClusters <-  sce_snp@colData %>% 
  as.data.table %>%
  .[, .(ID, cluster_hclust)]


fpcAssocs <- lapply(c(1:22, "X"), function(chr) {
  
  print(paste("chromosome",chr))
  
  results <- lapply(c(1:6), function(fpc) {
    
    # print(fpc)
    # slice <- 64
    file <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/results/chr",chr,"/chr",chr,"EUR.fpc",fpc,".glm.linear")
    
    fpcResults <- fread(file) %>%
      setnames(., "#CHROM", "CHR") %>%
      .[ID %in% pixClusters[,ID]] %>%
      .[, FPC := paste0("FPC",fpc)]
    
    return(fpcResults)
    
  }) %>%
    rbindlist
  
}) %>% 
  rbindlist


lapply(c(1:10), function(i){
  
  clustIDs <- pixClusters[cluster_hclust == i, ID]
  
  plotDT <- fpcAssocs[ID %in% clustIDs] %>%
    .[, association := sign(BETA)*log(P, 10)] %>%
    .[, sig := ifelse( P < 5E-8/6, "*", "")]
  
  plot <- ggplot(plotDT, aes(x = ID, y = FPC, fill = association)) +
    geom_tile() +
    geom_text(aes(label=sig))+
    scale_fill_gradient2() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90))
  
  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/plots/FPCassocs_pixelWiseSNPsCluster",i,".png"), width = 1200, height = 600)
  print(plot)
  dev.off()
  
  
  
})

medianFPCeffect <- lapply(c(1:10), function(i){
  
  clustIDs <- pixClusters[cluster_hclust == i, ID]
  
  clusEffect <- fpcAssocs[ID %in% clustIDs] %>%
    .[, .(medianBETA = median(BETA), medianP = median(P)), by = FPC] %>%
    .[, cluster := i]
  
  return(clusEffect)
  
}) %>%
  rbindlist  %>%
  .[, association := sign(medianBETA)*log(medianP, 10)]

plot <- ggplot(medianFPCeffect, aes(x = cluster, y = FPC, fill = association)) +
  geom_tile() +
  scale_fill_gradient2() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/plots/FPCassocsMedian_pixelWiseSNPsClusters.png", width = 1200, height = 600)
print(plot)
dev.off()


allFPCeffect <- lapply(c(1:10), function(i){
  
  clustIDs <- pixClusters[cluster_hclust == i, ID]
  
  clusEffect <- fpcAssocs[ID %in% clustIDs] %>%
    .[, cluster := i]
  
  return(clusEffect)
  
}) %>%
  rbindlist  %>%
  .[, log10P := (-1)*log(P, 10)] %>%
  .[, sig := ifelse(P < 5E-8/6, 1, 0)]

#    plot <- ggplot(allFPCeffect, aes(x = log10P, y = FPC)) +
#     geom_violin() +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 90)) +
#     facet_grid(cols = vars(cluster))

plot <- ggplot(allFPCeffect, aes(y = log10P, x = FPC)) +
  geom_violin() +
  geom_dotplot(binaxis='y', 
               stackdir='center', 
               dotsize = 1, aes(fill = as.factor(sig), alpha = 0.5)) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(cols = vars(cluster))

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/plots/FPCassocsAllSNPs_pixelWiseSNPsClusters.png", width = 1800, height = 1200)
print(plot)
dev.off()
