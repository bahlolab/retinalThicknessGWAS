## Enrichment analyses fPC reslts 
library(data.table)
library(magrittr)   
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)


annotResults <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/annotationsWithResults.csv")



## limit to genes in database
geneAnnot <- bitr(as.character(annotResults[,geneSYMBOL]), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

annotResultsDB <- merge.data.table(annotResults, geneAnnot, by.x = "geneSYMBOL", by.y= "SYMBOL") 


###########################################################################
## Run over-representation analyses.

## all fPC genes - use stricter top genes definition
fPCgenes <- lapply(c(1:6), function(i) {
  
  pCol <- paste0("P_FPC",i)
  print(paste(pCol))
  
  genes <- annotResultsDB[get(pCol) < 5E-8/6] %>%
    .[, .SD[which.max(evidenceScore)], by = ID] %>%
    .[, .(geneSYMBOL, ENTREZID)]
  
  
  return(genes)
}) %>%
  rbindlist



## run enrichment analysis  
## gene onology
go <- enrichGO(gene = unique(fPCgenes[,geneSYMBOL]),
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


## KEGG pathways
kegg <- enrichKEGG(gene = unique(fPCgenes[,ENTREZID]),
                   organism     = 'hsa',
                   keyType = 'ncbi-geneid',
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.1,
                   qvalueCutoff = 0.1)


if(nrow(kegg) > 0) { 
  
  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/overRepresentedKEGG_fPCall_topGenes.png"), width = 600, height = 600)
  p <- dotplot(kegg,showCategory = 20)  + 
    ggtitle("Over-represented KEGG") 
  print(p)
  dev.off()
  
}


## disease ontology
do <- enrichDO(gene = unique(fPCgenes[,ENTREZID]),
               pAdjustMethod = "BH",
               pvalueCutoff = 0.1,
               qvalueCutoff = 0.1)


if(nrow(do) > 0) {
  
  
  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/overRepresentedDO_fPCall_topGenes.png"), width = 600, height = 600)
  p <- dotplot(do,showCategory = 20)  + 
    ggtitle("Over-represented Disease Ontology") 
  print(p)
  dev.off()
  
}





## all Pixel-wise genes - top genes
pixelgenes <- annotResultsDB[P_allPixels < 5E-8/29041] %>%
  .[, .SD[which.max(evidenceScore)], by = ID] %>%
  .[, .(geneSYMBOL, ENTREZID)]



## run enrichment analysis  
## gene ontology
goPix <- enrichGO(gene = unique(pixelgenes[,geneSYMBOL]),
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


## KEGG pathwys
kegg <- enrichKEGG(gene = unique(pixelgenes[,ENTREZID]),
                   organism     = 'hsa',
                   keyType = 'ncbi-geneid',
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.1,
                   qvalueCutoff = 0.1)


if(nrow(kegg) > 0) { 
  
  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/overRepresentedKEGG_pixelsAll_topGenes.png"), width = 600, height = 600)
  p <- dotplot(kegg,showCategory = 20)  + 
    ggtitle("Over-represented KEGG") 
  print(p)
  dev.off()
  
}


## disease ontology
do <- enrichDO(gene = unique(pixelgenes[,ENTREZID]),
               pAdjustMethod = "BH",
               pvalueCutoff = 0.1,
               qvalueCutoff = 0.1)


if(nrow(do) > 0) {
  
  
  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/overRepresentedDO_pixelsAll_topGenes.png"), width = 600, height = 600)
  p <- dotplot(do,showCategory = 20)  + 
    ggtitle("Over-represented Disease Ontology") 
  print(p)
  dev.off()
  
}



## Now run comparison GO over-representation analysis,Gene Ontology only
symbolsCompare <- list(unique(pixelgenes[,geneSYMBOL]), unique(fPCgenes[,geneSYMBOL]))
names(symbolsCompare) <- c("pixelLevel", "fPC") 



compareGO <- compareCluster(symbolsCompare, 
                            fun="enrichGO",
                            ont = "ALL",
                            OrgDb = "org.Hs.eg.db",
                            keyType="SYMBOL",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.5)

dotplot(compareGO, showCategory = 15)


## and using less strict genes mapping
fPCgenesAll <- lapply(c(1:6), function(i) {
  
  pCol <- paste0("P_FPC",i)
  print(paste(pCol))
  
  genes <- annotResultsDB[get(pCol) < 5E-8/6] %>%
    .[, .(geneSYMBOL, ENTREZID)]
  
  
  return(genes)
}) %>%
  rbindlist

pixelgenesAll <- annotResultsDB[P_allPixels < 5E-8/29041] %>%
  .[, .(geneSYMBOL, ENTREZID)]

symbolsCompare <- list(unique(pixelgenesAll[,geneSYMBOL]), unique(fPCgenesAll[,geneSYMBOL]))
names(symbolsCompare) <- c("pixelLevel", "fPC") 



compareGO <- compareCluster(symbolsCompare, 
                            fun="enrichGO",
                            ont = "ALL",
                            OrgDb = "org.Hs.eg.db",
                            keyType="SYMBOL",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.5)

dotplot(compareGO, showCategory = 15)


## all together
symbolsCompare <- list( unique(pixelgenesAll[,geneSYMBOL]), unique(fPCgenesAll[,geneSYMBOL]),
                       unique(pixelgenes[,geneSYMBOL]), unique(fPCgenes[,geneSYMBOL]))
names(symbolsCompare) <- c("pixelLevel All Genes", "FPC All Genes",
                           "pixelLevel Top Genes", "FPC Top Genes") 



compareGO <- compareCluster(symbolsCompare, 
                            fun="enrichGO",
                            ont = "ALL",
                            OrgDb = "org.Hs.eg.db",
                            keyType="SYMBOL",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.5)

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/overRepresentedPathways_compareGO.png"), width = 900, height = 900)
dotplot(compareGO, showCategory = 15, label_format = 60)
ggtitle("Over-represented Pathways") %>%
  print
dev.off()

fwrite(as.data.table(compareGO), file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/overRepresentedPathways_compareGO.csv", sep = ",")
