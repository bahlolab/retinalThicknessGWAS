## Enrichment analyses fPC reslts 
library(data.table)
library(magrittr)   
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)


## For all idnetified sentinels, get results for all fPCs and top pixel.

pixelWiseSentinels <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/rsIDsSentinelsPixelwiseBonferroniSig.txt")

fpcSentinels <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/rsIDsSentinelsFPCsBonferroniSig.txt")

annot <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/geneSummaryAnnotations.csv")

fpcSNPs <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/FPCall/FUMA_job248242/snps.txt") %>%
  .[rsID==IndSigSNP] %>%
  .[, sentinel := rsID] %>%
  .[, .(sentinel, chr, pos)] %>%
  .[sentinel %in% fpcSentinels[,rsID]]

pixelSNPs <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/allPixels/FUMA_job248074/snps.txt") %>%
  .[rsID==IndSigSNP] %>%
  .[, sentinel := rsID] %>%
  .[, .(sentinel, chr, pos)] %>%
  .[sentinel %in% pixelWiseSentinels[,rsID]]

allSNPs <- rbind(fpcSNPs, pixelSNPs) %>%
  unique

annotPOS <- allSNPs[annot, on = "sentinel"] %>%
  .[, .(analysis = paste(sort(unique(analysis)), collapse = ";")), by = setdiff(names(.), "analysis")] %>%
  unique 


fpcResults <- lapply(c(1:22, "X"), function(x) {
  chrNo <- ifelse(x=="X", 23, x)
  
  print(paste("chr",x))
  
  chrPOS <- annotPOS[chr==chrNo]
  
  if(nrow(chrPOS) > 0) {
    
    fpcRes <- lapply(c(1:6), function(i) {
      print(paste(i))
      
      paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/results/chr",x,"/chr",x,"EUR.fpc",i,".glm.linear") %>%
        fread(., select = c("ID", "POS", "BETA", "P")) %>%
        setnames(., c("ID", "pos", "beta", "P")) %>%
        .[, FPC := paste0("FPC",i)] %>%
        .[pos %in% chrPOS[,pos]] %>%
        return(.)
    }) %>%
      rbindlist %>%
      .[, Padj := P * 6] %>%
      dcast(., ...  ~ FPC, value.var = c("beta", "P", "Padj")) %>%
      chrPOS[., on =  "pos"]
    
  } else{
    print(paste("no SNPs for chr",x))
    return(NULL)
  }
  
}) %>%
  rbindlist


pixelResultsPixelSentinels <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/allSentinelsAllPixelsResults_clumpThresh0.001_withOverlap.csv") %>%
  .[, Padj := P*29041] %>%
  .[Padj < 5E-8] %>%
  .[, .SD[which.min(Padj)], by = ID] %>%
  .[, .(chrom, POS, ID, BETA, P, Padj)] %>%
  setnames(., c("chr", "pos", "ID", "beta_allPixels", "P_allPixels", "Padj_allPixels"))

pixels <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/pixels.txt") %>%
  setnames(., c("pixel", "y", "x"))

pixelResultsFPCsentinels <- lapply(c(1:22, "X"), function(chr) {
  
  ## full sig results
  dir <- paste0("/vast/scratch/users/jackson.v/retThickness/GWASfollowup/pixelWiseResultsFPCsentinels/results/chr",chr)
  print(paste("chr",chr))
  
  chrResults <- lapply(c(1:119), function(slice) {
    
    chrPOS <- annotPOS[chr==chr]
    
    file <- paste0(dir,"/chr",chr,"Slice",slice,"_sentinels.txt")
    
    sliceResults <- fread(file) %>%
      setnames(., "#POS", "POS") %>%
      .[POS %in% chrPOS[,pos]] %>%
      .[, chrom := chr] %>%
      .[, Padj := P*29041] %>%
      .[, .(chrom, POS, ID, BETA, P, Padj)] %>%
      setnames(., c("chr", "pos", "ID", "beta_allPixels", "P_allPixels", "Padj_allPixels"))
    
    return(sliceResults)
    
  }) %>%
    rbindlist %>%
    .[, .SD[which.min(Padj_allPixels)], by = ID] 
  
  return(chrResults)
  
}) %>%
  rbindlist

pixelResults <- rbind(pixelResultsPixelSentinels, pixelResultsFPCsentinels) %>%
  unique %>%
  .[, chr := ifelse(chr=="X", 23, as.integer(chr))]

annotResults <- fpcResults %>%
  pixelResults[., on = c("ID", "chr", "pos")] 




## limit to genes in database
geneAnnot <- bitr(as.character(annotResults[,gene]), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

annotResultsDB <- merge.data.table(annotResults, geneAnnot, by.x = "gene", by.y= "SYMBOL") %>%
  .[, mappingSum := sum_cols-(OMIM + mousePhenotype)]



###########################################################################
## Run over-representation analyses.

## all fPC genes - use stricter top genes definition
fPCgenes <- lapply(c(1:6), function(i) {
  
  pCol <- paste0("P_FPC",i)
  print(paste(pCol))
  
  genes <- annotResultsDB[get(pCol) < 5E-8/6] %>%
    .[, .SD[which.max(mappingSum)], by = sentinel] %>%
    .[, .(gene, ENTREZID)]
  
  
  return(genes)
}) %>%
  rbindlist



## run enrichment analysis  
## gene onology
go <- enrichGO(gene = unique(fPCgenes[,gene]),
               ont = "ALL",
               OrgDb = "org.Hs.eg.db",
               keyType="SYMBOL",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.5)




png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/overRepresentedPathways_fPCall_topGenes.png"), width = 600, height = 600)
dotplot(go,showCategory = 20)  + 
  ggtitle("Over-represented Pathways") %>%
  print
dev.off()

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/overRepresentedPathwaysGenes_fPCall_topGenes.png"), width = 600, height = 600)
heatplot(go,showCategory = 20) + 
  ggtitle("Genes in Over-represented Pathways") %>%
  print
dev.off()

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/overRepresentedPathwaysRelationships_fPCall_topGenes.png"), width = 600, height = 600)
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
  
  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/overRepresentedKEGG_fPCall_topGenes.png"), width = 600, height = 600)
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
  
  
  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/overRepresentedDO_fPCall_topGenes.png"), width = 600, height = 600)
  p <- dotplot(do,showCategory = 20)  + 
    ggtitle("Over-represented Disease Ontology") 
  print(p)
  dev.off()
  
}



## ORA gene ontology for each fPC separately.

goResults <- lapply(c(1:6), function(i) {
  
  
  pCol <- paste0("P_FPC",i)
  print(paste(pCol))
  
  fPCgenes <- fpcResults[get(pCol) < 5E-8/6, gene]
  
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
  
  
  
  
  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/overRepresentedPathways_fPC",i,".png"), width = 600, height = 600)
  dotplot(go,showCategory = 20)  + 
    ggtitle("Over-represented Pathways") %>%
    print
  dev.off()
  
  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/overRepresentedPathwaysGenes_fPC",i,".png"), width = 600, height = 600)
  heatplot(go,showCategory = 20) + 
    ggtitle("Genes in Over-represented Pathways") %>%
    print
  dev.off()
  
  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/overRepresentedPathwaysRelationships_fPC",i,".png"), width = 600, height = 600)
  emapplot(pairwise_termsim(go),showCategory = 20, layout='nicely', cex_label_category = 0.7) + 
    ggtitle("Relationship between Over-represented Pathways") %>%
    print
  dev.off()
  
  return(go)
})





## all Pixel-wise genes - top genes
pixelgenes <- annotResultsDB[P_allPixels < 5E-8/29041] %>%
  .[, .SD[which.max(mappingSum)], by = sentinel] %>%
  .[, .(gene, ENTREZID)]



## run enrichment analysis  
## gene ontology
goPix <- enrichGO(gene = unique(pixelgenes[,gene]),
                  ont = "ALL",
                  OrgDb = "org.Hs.eg.db",
                  keyType="SYMBOL",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.5)




png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/overRepresentedPathways_pixelsAll_topGenes.png"), width = 600, height = 600)
dotplot(goPix,showCategory = 20)  + 
  ggtitle("Over-represented Pathways") %>%
  print
dev.off()

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/overRepresentedPathwaysGenes_pixelsAll_topGenes.png"), width = 600, height = 600)
heatplot(goPix,showCategory = 20) + 
  ggtitle("Genes in Over-represented Pathways") %>%
  print
dev.off()

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/overRepresentedPathwaysRelationships_pixelsAll_topGenes.png"), width = 600, height = 600)
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
  
  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/overRepresentedKEGG_pixelsAll_topGenes.png"), width = 600, height = 600)
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
  
  
  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/overRepresentedDO_pixelsAll_topGenes.png"), width = 600, height = 600)
  p <- dotplot(do,showCategory = 20)  + 
    ggtitle("Over-represented Disease Ontology") 
  print(p)
  dev.off()
  
}



## Now run comparison GO over-representation analysis,Gene Ontology only
symbolsCompare <- list(unique(pixelgenes[,gene]), unique(fPCgenes[,gene]))
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
    .[, .(gene, ENTREZID)]
  
  
  return(genes)
}) %>%
  rbindlist

pixelgenesAll <- annotResultsDB[P_allPixels < 5E-8/29041] %>%
  .[, .(gene, ENTREZID)]

symbolsCompare <- list(unique(pixelgenesAll[,gene]), unique(fPCgenesAll[,gene]))
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
