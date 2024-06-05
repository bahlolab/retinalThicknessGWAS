library(data.table)
library(magrittr)
library(funData)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(UpSetR)
library(stringr)
library(DescTools)
library(patchwork)
library(ggside)
require(gridExtra)
library(purrr)


# chr <- 9

# sentinels <- lapply(c(1:22), function(chr) {

# chrSent <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/chr",chr,"sentinelsIDonly_clumpThresh0.001_withOverlap.txt") %>%
#     fread(., header=F)
# return(chrSent)
# }) %>%
# rbindlist %>%
# setnames(., "ID")

# ## read in results for all sentinels
# results <- lapply(c(1:22, "X"), function(chr) {

#   dir <- paste0("/vast/scratch/users/jackson.v/retThickness/GWAS/GWsigResults/chr",chr)

# print(paste("chromosome",chr))
# chrResults <- lapply(c(1:119), function(slice) {

#   print(paste(slice))

#   file <- paste0(dir,"/chr",chr,"Slice",slice,"_5e-5Sig.txt")

#   sliceResults <- fread(file) %>%
#     setnames(., "#POS", "POS") %>%
#     .[ID %in% sentinels[,ID]]
#   return(sliceResults)

# }) %>%
# rbindlist %>%
# .[, chrom := chr]

# return(chrResults)
# }) %>%
# rbindlist   

# fwrite(results, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/allSentinelsAllPixelsResults_clumpThresh0.001_withOverlap.csv", sep = ",")

### Plots showing distribution of loci, by pixel
results <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/allSentinelsAllPixelsResults_clumpThresh0.001_withOverlap.csv")
load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/rawData/example_data.RData")


dt <- results[P<5e-8/29041] %>%
  .[ , .N, by = pixel] %>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert = T)] %>%
  .[!is.na(y)]

plot <- ggplot(dt) +
  geom_tile(aes(x = x, y = y, fill = N)) +
  scale_fill_gradient2() +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 16))+
  geom_path(aes(x=col,y=diag1),color="#4B4B4B", data=res[res$variable=="logFC",])+
  geom_path(aes(x=col,y=diag2),color="#4B4B4B", data=res[res$variable=="logFC",])+
  geom_path(aes(x=col,y=line1),color="#4B4B4B", data=res[res$variable=="logFC",])+
  geom_path(aes(x=col,y=line2),color="#4B4B4B", data=res[res$variable=="logFC",])+
  geom_path(aes(x=col,y=row,group=area),color="#4B4B4B",data = areas[areas$edtrs,],size = 0.5)   

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/plots/allChr_pixelsPerSentinel_clumpThresh0.001_BonfSig_grid.png"), width = 1200, height = 1200)
print(plot)
dev.off()


png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/plots/allChr_pixelsPerSentinel_clumpThresh0.001_BonfSig_hist.png"), width = 800, height = 800)
ggplot(dt, aes(x=N)) +
  geom_histogram()
dev.off()


## Upsettr plot for fPCs only.
fpcSentinels <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/sentinels/allChr_sentinel_clumpThresh0.001_withOverlap.txt") %>%
  .[nSNPsLocus >= 5]

loci <- lapply(c(1:22, "X"), function(chr) {
  
  print(paste("chromosome",chr))
  
  results <- lapply(c(1:6), function(fpc) {
    
    # print(fpc)
    # slice <- 64
    file <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/results/chr",chr,"/chr",chr,"EUR.fpc",fpc,".glm.linear")
    
    fpcResults <- fread(file) %>%
      setnames(., "#CHROM", "CHR") %>%
      .[ID %in% fpcSentinels[,ID]] %>%
      .[, FPC := paste0("FPC",fpc)]
    
    return(fpcResults)
    
  }) %>%
    rbindlist
  
}) %>% rbindlist

lociUpsett <- loci[, FPCsig := as.integer(P < 5E-8/6)] %>%
  dcast(., ID ~ FPC, value.var = "FPCsig")

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/plots/upsetFPCs.png", width = 1200, height = 600)
upset(lociUpsett, nsets = 6, text.scale = 2) %>% print
dev.off()



loci[P<5e-8/6] %>%
  .[ , .N, by = FPC] 




## Summary plot for main paper 
annot <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/geneSummaryAnnotations.csv")

pixelWiseSentinelsFull <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/cojoAnalysis/output/bonfSigLociWithSecondary.csv")
sents <- pixelWiseSentinelsFull[,.(CHR, POS, ID, rsID)]
secondary <- pixelWiseSentinelsFull[,.(CHR, POS_secondary, SNP_secondary, SNP_secondary)] %>%
  na.omit %>%
  setnames(., names(sents))
pixelWiseSentinels <- rbind(sents, secondary)

fpcSentinelsFull <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/cojoAnalysis/output/bonfSigFPCLociWithSecondary.csv")
sents <- fpcSentinelsFull[,.(CHR, POS, ID, rsID)]
secondary <- fpcSentinelsFull[,.(CHR, POS_secondary, SNP_secondary, SNP_secondary)] %>%
  na.omit %>%
  setnames(., names(sents))
fpcSentinels <- rbind(sents, secondary)




fpcSNPs <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/FPCall/FUMA_job483186/snps.txt") %>%
  .[rsID==IndSigSNP] %>%
  .[, sentinel := rsID] %>%
  .[, .(sentinel, chr, pos)] %>%
  .[sentinel %in% fpcSentinels[,rsID]]

pixelSNPs <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/allPixels/FUMA_job483192/snps.txt") %>%
  .[rsID==IndSigSNP] %>%
  .[, sentinel := rsID] %>%
  .[, .(sentinel, chr, pos)] %>%
  .[sentinel %in% pixelWiseSentinels[,rsID]]

allSNPs <- rbind(pixelWiseSentinels, fpcSentinels) %>%
  unique

fwrite(allSNPs[,.(ID)], file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/misc/allSigSNPs.txt")


annotPOS <- allSNPs[annot, on = c("ID" = "sentinel")] %>%
  .[, .(analysis = paste(sort(unique(analysis)), collapse = ";")), by = setdiff(names(.), "analysis")] %>%
  unique %>%
  setnames(., c("CHR", "POS"), c("chr", "pos"))


fpcResults <- lapply(c(1:22, "X"), function(x) {
  chrNo <- ifelse(x=="X", 23, x)
  
  print(paste("chr",x))
  
  chrPOS <- annotPOS[chr==chrNo]
  
  if(nrow(chrPOS) > 0) {
    
    fpcRes <- lapply(c(1:6), function(i) {
      print(paste(i))
      
      paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/results/chr",x,"/chr",x,"EUR.fpc",i,".glm.linear") %>%
        fread(., select = c("ID", "POS", "A1", "BETA", "P")) %>%
        setnames(., c("ID", "pos", "A1", "beta", "P")) %>%
        .[, FPC := paste0("FPC",i)] %>%
        .[pos %in% chrPOS[,pos]] %>%
        return(.)
    }) %>%
      rbindlist %>%
      .[, Padj := P*6] %>%
      dcast(., ...  ~ FPC, value.var = c("beta", "P", "Padj")) %>%
      chrPOS[., on =  c("ID", "pos")]
    
  } else{
    print(paste("no SNPs for chr",x))
    return(NULL)
  }
  
}) %>%
  rbindlist %>%
  .[, chr := ifelse(chr=="X", 23, as.integer(chr))]


# pixelResultsPixelSentinels <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/allSentinelsAllPixelsResults_clumpThresh0.001_withOverlap.csv") %>%
#   .[, Padj := P*29041] %>%
#   .[Padj < 5E-8] %>%
#   .[, .SD[which.min(Padj)], by = ID] %>%
#   .[, .(chrom, POS, A1, ID, BETA, P, Padj)] %>%
#   setnames(., c("chr", "pos", "A1", "ID", "beta_allPixels", "P_allPixels", "Padj_allPixels"))
# 
# pixels <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/pixels.txt") %>%
#   setnames(., c("pixel", "y", "x"))
# 
# pixelResultsFPCsentinels <- lapply(c(1:22, "X"), function(chr) {
#   
#   ## full sig results
#   dir <- paste0("/vast/scratch/users/jackson.v/retThickness/GWASfollowup/pixelWiseResultsFPCsentinels/results/chr",chr)
#   print(paste("chr",chr))
#   
#   chrResults <- lapply(c(1:119), function(slice) {
#     
#     chrPOS <- annotPOS[chr==chr]
#     
#     file <- paste0(dir,"/chr",chr,"Slice",slice,"_sentinels.txt")
#     
#     sliceResults <- fread(file) %>%
#       setnames(., "#POS", "POS") %>%
#       .[POS %in% chrPOS[,pos]] %>%
#       .[, chrom := chr] %>%
#       .[, Padj := P*29041] %>%
#       .[, .(chrom, POS, A1, ID, BETA, P, Padj)] %>%
#       setnames(., c("chr", "pos", "A1", "ID", "beta_allPixels", "P_allPixels", "Padj_allPixels"))
#     
#     return(sliceResults)
#     
#   }) %>%
#     rbindlist %>%
#     .[, .SD[which.min(Padj_allPixels)], by = ID] 
#   
#   return(chrResults)
#   
# }) %>%
#   rbindlist


pixelResults <- lapply(c(1:22,  "X"), function(chr) {
  
  dir <- paste0("/vast/scratch/users/jackson.v/retThickness/GWASfollowup/pixelWiseResultsAllSentinels/results/chr",chr)
  print(paste("chr",chr))
  
  chrResults <- lapply(c(1:119), function(slice) {
    
    
    file <- paste0(dir,"/chr",chr,"Slice",slice,"_sentinels.txt")
    
      sliceResults <- fread(file) %>%
        setnames(., "#POS", "POS") %>%
        .[, chrom := chr] %>%
        .[, Padj := as.numeric(P)*29041] %>%
        .[, .(chrom, POS, A1, ID, BETA, P, Padj)] %>%
        setnames(., c("chr", "pos", "A1", "ID", "beta_allPixels", "P_allPixels", "Padj_allPixels"))
        
      return(sliceResults)
    }) %>%
    rbindlist %>%
     .[, .SD[which.min(Padj_allPixels)], by = ID] 
   
    return(chrResults) 
  }) %>%
    rbindlist  %>%
  unique %>%
  .[, chr := ifelse(chr=="X", 23, as.integer(chr))] %>%
  .[, pos := as.integer(pos)] %>%
  .[, beta_allPixels := as.numeric(beta_allPixels)] %>%
  .[, P_allPixels := as.numeric(P_allPixels)]

annotResults <- fpcResults %>%
  pixelResults[., on = c("ID", "chr", "pos")] %>%
  .[, gene := paste0(gene," (",rsID,")")]


## for plot, select gene(s) with most lines fof evidence for each sentinel.
annotPlot <- annotResults %>%
  .[, .SD[which.max(evidenceScore)], by = ID] 





## extract results for each sentinel
analysesCols <- names(annotPlot)[names(annotPlot) %like any% c("Padj_%", "beta_%")]
betaCols <- names(annotPlot)[names(annotPlot) %like% "beta_%"]
PadjCols <- names(annotPlot)[names(annotPlot) %like% "Padj_%"]

betasLong <- annotPlot %>%
  reshape2::melt(., id.vars = c("rsID", "gene"), 
                 measure.vars = betaCols,
                 variable.name = "analysis", value.name = "beta") %>%
  as.data.table %>%
  .[, analysis := str_split(analysis, "_", simplify = TRUE)[, 2]]

pAdjLong <- annotPlot %>%
  reshape2::melt(., id.vars = c("rsID", "gene"), 
                 measure.vars = PadjCols,
                 variable.name = "analysis", value.name = "Padj") %>%
  as.data.table %>%
  .[, analysis := str_split(analysis, "_", simplify = TRUE)[, 2]]


## put all results together
plotAnalyses <-  betasLong[pAdjLong, on = c("rsID", "gene", "analysis")] %>%
  .[, association := sign(beta)*log(Padj, 10)]  %>%
  .[, .(gene, analysis, beta, Padj, association)] %>%
  .[, cat:= ifelse(analysis=="allPixels", "pixelwise", "FPC")] %>%
  na.omit


## order by adjusted p-value
geneOrder <- plotAnalyses %>%
  .[, .SD[which.min(Padj)], by = c("gene")] %>%
  .[order(Padj)] %>%
  .[,gene] 
## dt for plotting annotations
geneAnnotPlot <- melt(annotResults, id.vars = c("rsID", "gene", "analysis"), measure.vars = c("proximity", "exonic", "CADD20", "eQTLBlood", "eQTLBrain", "eQTLRetina", "chromatinInteractionCortex", "chromatinInteractionRetina", "OMIM", "mousePhenotype")) %>%
  .[, evidence := ifelse(value ==1, "Y", "N") %>% as.factor] 
annotOrder <- c("proximity", "exonic", "CADD20", "eQTLBlood", "eQTLBrain", "eQTLRetina", "chromatinInteractionCortex","chromatinInteractionRetina", "OMIM", "mousePhenotype")


## plot - split into 
colours <- rev(brewer.pal(9,"RdBu"))
biColours <- brewer.pal(3,"Purples")[c(1,3)]
lim <- plotAnalyses[,association] %>% abs %>% max(., na.rm=T) %>% ceiling()

to <- seq(0, length(geneOrder), length.out = 5) %>% round
from <- shift(to) + 1
# combine the "from" and "to" sequences into a data.table
index <- data.table(from = from[-1], to = to[-1])

plots <- lapply(c(1:4), function(i) {
  
  fr <- index[i,from]
  to <- index[i, to]
  
  plotAnalyses %>% 
    ggplot(., aes(x = analysis, y = gene)) +
    geom_tile(aes(fill = association), color = "white", lwd = 0.3,linetype = 1) +
    #  scale_fill_manual(values = color_scale) +
    labs(y = "Genes", title = "", x = "") +
    scale_fill_gradientn(colors=colours, limits = c(-lim, lim), values = c(0, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 1) ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_y_discrete(limits = rev(geneOrder[fr:to]) ) +
    theme(legend.position = "top", legend.key.width= unit(0.8, 'cm')) +
    geom_ysidetile(data = geneAnnotPlot, aes(x = variable, yfill = evidence), color = "white", lwd = 0.3,linetype = 1) +
    theme(ggside.panel.scale = 1) +
    scale_yfill_manual(values = biColours) +
    scale_ysidex_discrete(limits = annotOrder)  %>%
    return
  
})

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/prioritisedGenes_20240605.png", width = 1800, height = 1100)
reduce(plots, `|`) + plot_layout(guides = "collect") & theme(legend.position = "top") 
dev.off()

assocOrder <- c("allPixels", paste0("FPC", 1:6))

plots <- lapply(c(1:4), function(i) {
  
  fr <- index[i,from]
  to <- index[i, to]
  
  plotAnalyses %>% 
    ggplot(., aes(y = analysis, x = gene)) +
    geom_tile(aes(fill = association), color = "white", lwd = 0.3,linetype = 1) +
    #  scale_fill_manual(values = color_scale) +
    labs(x = "Genes", title = "", y = "") +
    scale_fill_gradientn(colors=colours, limits = c(-lim, lim), values = c(0, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 1) ) +
    theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=1)) +
    scale_x_discrete(limits = geneOrder[fr:to] , position = "top") +
    scale_y_discrete(limits = rev(assocOrder)) +
    theme(legend.position = "bottom", legend.key.width= unit(0.8, 'cm')) +
    ggside(x.pos = "bottom") +
    geom_xsidetile(data = geneAnnotPlot, aes(y = variable, xfill = evidence), color = "white", lwd = 0.3,linetype = 1) +
    theme(ggside.panel.scale = 1) +
    scale_xfill_manual(values = biColours) +
    scale_xsidey_discrete(limits = rev(annotOrder))  %>%
    return
  
})

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/prioritisedGenes_portrait_20240605.png", width = 1600, height = 1600)
reduce(plots, `/`) + plot_layout(guides = "collect") & theme(legend.position = "bottom") 
dev.off()








## Enrichment analyses fPC reslts 
library(data.table)
library(magrittr)   
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)


## over-representation analyses






## all fPC genes
fPCgenes <- lapply(c(1:6), function(i) {
  
  pCol <- paste0("P_FPC",i)
  print(paste(pCol))
  
  genes <- annotResults[get(pCol) < 5E-8/6, .(gene)] %>%
    .[, c("geneSYMBOL", "snp") := tstrsplit(gene, " ")] %>%
    .[,geneSYMBOL]
  
  return(genes)
}) %>%
  unlist %>%
  unique
  

## limit to genes in database
orgdb <- org.Hs.eg.db # change to human DBs
cols <- c( "ENTREZID")

geneAnnotExtra <- AnnotationDbi::select(orgdb, keys=as.character(fPCgenes), columns=cols, keytype="SYMBOL") %>%
  as.data.table %>%
  na.omit


## run enrichment analysis  
go <- enrichGO(gene = unique(geneAnnotExtra[,SYMBOL]),
               ont = "ALL",
               OrgDb = "org.Hs.eg.db",
               keyType="SYMBOL",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.5)




png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/overRepresentedPathways_fPCall.png"), width = 600, height = 600)
dotplot(go,showCategory = 20)  + 
  ggtitle("Over-represented Pathways") %>%
  print
dev.off()

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/overRepresentedPathwaysGenes_fPCall.png"), width = 600, height = 600)
heatplot(go,showCategory = 20) + 
  ggtitle("Genes in Over-represented Pathways") %>%
  print
dev.off()

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/overRepresentedPathwaysRelationships_fPCall.png"), width = 600, height = 600)
emapplot(pairwise_termsim(go),showCategory = 20, layout='nicely', cex_label_category = 0.7) + 
  ggtitle("Relationship between Over-represented Pathways") %>%
  print
dev.off()


## all Pixel-wise genes
pixelgenes <- annotResults[P_allPixels < 5E-8/29041, .(gene)] %>%
  .[, c("geneSYMBOL", "snp") := tstrsplit(gene, " ")] %>%
  .[,geneSYMBOL] %>%
  unique


## limit to genes in database
orgdb <- org.Hs.eg.db # change to human DBs
cols <- c( "ENTREZID")

geneAnnotExtra <- AnnotationDbi::select(orgdb, keys=as.character(pixelgenes), columns=cols, keytype="SYMBOL") %>%
  as.data.table %>%
  na.omit


## run enrichment analysis  
goPix <- enrichGO(gene = unique(geneAnnotExtra[,SYMBOL]),
                  ont = "ALL",
                  OrgDb = "org.Hs.eg.db",
                  keyType="SYMBOL",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.5)




png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/overRepresentedPathways_pixelsAll.png"), width = 600, height = 600)
dotplot(goPix,showCategory = 20)  + 
  ggtitle("Over-represented Pathways") %>%
  print
dev.off()

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/overRepresentedPathwaysGenes_pixelsAll.png"), width = 600, height = 600)
heatplot(goPix,showCategory = 20) + 
  ggtitle("Genes in Over-represented Pathways") %>%
  print
dev.off()

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/overRepresentedPathwaysRelationships_pixelsAll.png"), width = 600, height = 600)
emapplot(pairwise_termsim(goPix),showCategory = 20, layout='nicely', cex_label_category = 0.7) + 
  ggtitle("Relationship between Over-represented Pathways") %>%
  print
dev.off()











## limit to genes in database
annotResults <- annotResults %>%
  .[, c("gene", "snp") := tstrsplit(gene, " ")]
geneAnnot <- bitr(as.character(annotResults[,gene]), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

annotResultsDB <- merge.data.table(annotResults, geneAnnot, by.x = "gene", by.y= "SYMBOL") 

## all fPC genes stricter top genes
fPCgenes <- lapply(c(1:6), function(i) {
  
  pCol <- paste0("P_FPC",i)
  print(paste(pCol))
  
  genes <- annotResultsDB[get(pCol) < 5E-8/6] %>%
    .[, .SD[which.max(evidenceScore)], by = ID] %>%
    .[, .(gene, ENTREZID)]
  
  
  return(genes)
}) %>%
  rbindlist %>%
  unique



## run enrichment analysis  
go <- enrichGO(gene = unique(fPCgenes[,gene]),
               ont = "ALL",
               OrgDb = "org.Hs.eg.db",
               keyType="SYMBOL",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.5)




png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/overRepresentedPathways_fPCall_topGenes.png"), width = 600, height = 600)
dotplot(go,showCategory = 20)  + 
  ggtitle("Over-represented Pathways") %>%
  print
dev.off()

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/overRepresentedPathwaysGenes_fPCall_topGenes.png"), width = 600, height = 600)
heatplot(go,showCategory = 20) + 
  ggtitle("Genes in Over-represented Pathways") %>%
  print
dev.off()

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/overRepresentedPathwaysRelationships_fPCall_topGenes.png"), width = 600, height = 600)
emapplot(pairwise_termsim(go),showCategory = 20, layout='nicely', cex_label_category = 0.7) + 
  ggtitle("Relationship between Over-represented Pathways") %>%
  print
dev.off()







## all Pixel-wise genes top genes
pixelgenes <- annotResultsDB[P_allPixels < 5E-8/29041] %>%
  .[, .SD[which.max(evidenceScore)], by = ID] %>%
  .[, .(gene, ENTREZID)] %>%
  unique




## run enrichment analysis  
goPix <- enrichGO(gene = unique(pixelgenes[,gene]),
                  ont = "ALL",
                  OrgDb = "org.Hs.eg.db",
                  keyType="SYMBOL",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.5)




png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/overRepresentedPathways_pixelsAll_topGenes.png"), width = 600, height = 600)
dotplot(goPix,showCategory = 20)  + 
  ggtitle("Over-represented Pathways") %>%
  print
dev.off()

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/overRepresentedPathwaysGenes_pixelsAll_topGenes.png"), width = 600, height = 600)
heatplot(goPix,showCategory = 20) + 
  ggtitle("Genes in Over-represented Pathways") %>%
  print
dev.off()

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/overRepresentedPathwaysRelationships_pixelsAll_topGenes.png"), width = 600, height = 600)
emapplot(pairwise_termsim(goPix),showCategory = 20, layout='nicely', cex_label_category = 0.7) + 
  ggtitle("Relationship between Over-represented Pathways") %>%
  print
dev.off()






























symbolsCompare <- list(unique(pixelgenes[,gene]), unique(fPCgenes[,gene]))
names(symbolsCompare) <- c("pixelLevel", "fPC") 


entrezCompare <- list(unique(pixelgenes[,ENTREZID]), unique(fPCgenes[,ENTREZID]))
names(entrezCompare) <- c("pixelLevel", "fPC") 

compareGO <- compareCluster(symbolsCompare, 
                            fun="enrichGO",
                            ont = "ALL",
                            OrgDb = "org.Hs.eg.db",
                            keyType="SYMBOL",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.5)

dotplot(compareGO, showCategory = 15)


compareKEGG <- compareCluster(entrezCompare, 
                              fun="enrichKEGG",
                              organism  = 'hsa',
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.5,
                              qvalueCutoff = 0.5)

dotplot(compareKEGG, showCategory = 20)



compareDO <- compareCluster(entrezCompare, 
                            fun="enrichDO",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.5,
                            qvalueCutoff = 0.5)

dotplot(compareDO, showCategory = 10)



compareReactome <- compareCluster(entrezCompare, 
                                  fun="enrichPathway",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff = 0.5,
                                  qvalueCutoff = 0.5)

dotplot(compareReactome, showCategory = 10)






dgnvFPC <- enrichDGNv(fpcSentinels[,rsID],
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.5,
                      qvalueCutoff = 0.5)

dgnvPixel <- enrichDGNv(pixelWiseSentinels[,rsID],
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.5,
                        qvalueCutoff = 0.5)

dotplot(dgnvFPC, showCategory = 20)
dotplot(dgnvPixel, showCategory = 20)




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
