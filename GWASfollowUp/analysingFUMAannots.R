#!/usr/bin/env Rscript


library(data.table)
library(magrittr)
library(tidyverse)
library(optparse)
library(DescTools)
library(stringr)
library(RColorBrewer)
# option_list <-  list(
#   make_option(c("-e", "--eqtl"), type="character", default=NULL,
#               help="eQTL file", metavar="character"),
#   make_option(c("-g", "--genes"), type="character", default=NULL,
#               help="eQTL file", metavar="character"),
#   make_option(c("-s", "--snps"), type="character", default=NULL,
#               help="snps file", metavar="character"),
#   make_option(c("-cis", "--cisnps"), type="character", default=NULL,
#               help="ciSNPs file", metavar="character"),
#   make_option(c("-ci", "--ci"), type="character", default=NULL,
#               help="ci file", metavar="character"),
#   make_option(c("-gc", "--gwascat"), type="character", default=NULL,
#               help="gwas catalogue file", metavar="character"),
# );
# 
# opt_parser <- OptionParser(option_list=option_list)
# opt <- parse_args(opt_parser)
# 
# eqtlFile <-  opt$eqtl
# genesFile <- opt$genes
# snpsFile <- opt$snps
# ciSNPsFile <- opt$cisnps
# ciFile <- opt$ci
# gwascatFile <- opt$gwascat
# 
# eqtlFile <-  "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/eqtl.txt"
# genesFile <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/genes.txt"
# snpsFile <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/snps.txt"
# ciSNPsFile <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/ciSNPs.txt"
# ciFile <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/ci.txt"
# gwascatFile <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/gwascatalog.txt"


# annot <- lapply(c("pixelWise", "FPC1", "FPC2", "FPC3", "FPC4", "FPC5", "FPC6"), function(analysis) {
annot <- lapply(c("pixelWise", "FPCall"), function(analysis) {
    
  if(analysis=="pixelWise") { 
    dir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/allPixels/FUMA_job248074"}

  if(analysis=="FPCall") { 
    dir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/FPCall/FUMA_job248242"}
  
  if(analysis=="FPC1") { 
    dir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/FPC1/FUMA_job246508"}
  
  if(analysis=="FPC2") { 
dir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/FPC2/FUMA_job246510"}
  
  if(analysis=="FPC3") {
dir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/FPC3"}
  
  if(analysis=="FPC4") {
dir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/FPC4/FUMA_job246513"}
  
  if(analysis=="FPC5") {
dir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/FPC5/FUMA_job246514"}
  
  if(analysis=="FPC6") {
dir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/FPC6/FUMA_job246516"}

eqtlFile <-  paste0(dir,"/eqtl.txt")
genesFile <- paste0(dir,"/genes.txt")
snpsFile <- paste0(dir,"/snps.txt")
ciSNPsFile <- paste0(dir,"/ciSNPs.txt")
ciFile <- paste0(dir,"/ci.txt")
gwascatFile <- paste0(dir,"/gwascatalog.txt")

eqtl <- fread(eqtlFile)
genes <- fread(genesFile)
snps <- fread(snpsFile)

retinaEQTLs <- snps[, .(IndSigSNP, uniqID, r2)] %>%
    .[eqtl[db == "EyeGEx", .(uniqID, symbol, p, FDR)], on = "uniqID"] %>%
    .[, .SD[which.min(p)], by = c("symbol", "IndSigSNP")]

allEQTLs <- snps[, .(IndSigSNP, uniqID, r2)] %>%
    .[eqtl[, .(uniqID, tissue, symbol, p, FDR)], on = "uniqID"] %>%
    .[, .SD[which.min(p)], by = c("symbol", "IndSigSNP", "tissue")] %>%
    .[r2 > 0.5] %>%
    setnames(., c("gene", "sentinel", "eQTLtissue", "mappedSNP", "r2", "P", "FDR")) %>%
    .[gene %in% genes[,symbol]]

fwrite(allEQTLs, file = paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/",analysis,"_eQTLs.csv"), sep = ",")

ciSNPs <- fread(ciSNPsFile) %>%
.[snps[, .(IndSigSNP, uniqID, r2)] , on = "uniqID"] %>%
.[!is.na(rsID)] %>%
.[, .(rsID, IndSigSNP, uniqID, r2)] %>%
unique

ci <- fread(ciFile) %>%
    .[DB == "Giusti-Rodriguez_et_al_2019"] %>%
    .[!is.na(genes)] %>%
    .[, .(rsID = unlist(strsplit(SNPs, ";"))), by = c("region1", "region2", "FDR",  "tissue/cell", "genes")] %>%  
    .[, .(gene = unlist(strsplit(genes, ":"))), by = c("region1", "region2", "FDR",  "tissue/cell", "rsID")] %>%
    .[ciSNPs, on = "rsID"] %>%
    .[r2 > 0.5 ] %>%
    .[!is.na(FDR)] %>%
    .[genes[,.(ensg, symbol)], on = c("gene" = "ensg")] %>%
    .[!is.na(IndSigSNP)] %>%    
    .[, c("IndSigSNP", "rsID", "r2", "FDR", "tissue/cell", "symbol")] %>%
    setnames(., "tissue/cell" , "tissue") %>%
    setnames(., "symbol" , "gene") %>%
     setnames(., "rsID", "mappedSNP") %>%
    setnames("IndSigSNP", "sentinel") %>%
    .[gene %in% genes[,symbol]]

fwrite(ci, file = paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/",analysis,"_chromatinInteractions.csv"), sep = ",")


position <- snps[r2 > 0.5] %>%
    .[, .(uniqID, rsID, r2, IndSigSNP, nearestGene, dist, func, CADD)]



sents <- position[,IndSigSNP] %>% unique

sentsAnnot <- lapply(sents, function(s) {

     # print(paste(s))

    sentsMaps <- position[IndSigSNP == s]

    inGene <- sentsMaps %>%
        tidyr::separate_rows(., nearestGene, dist, sep = ":") %>%
        as.data.table() %>%
        .[dist==0] %>%
        .[, .SD[which.max(r2)], by = c("IndSigSNP", "nearestGene")] %>%
        .[nearestGene %in% genes[,symbol]] %>%
        .[, exonic := ifelse(func == "exonic", 1, 0 )] %>%
        .[, CADD20 := ifelse(CADD >= 20, 1, 0)]


    if(nrow(inGene) > 0) {
        return(inGene)
    } else {
        sentGene <- sentsMaps %>%
            tidyr::separate_rows(., nearestGene, dist, sep = ":") %>%
            as.data.table() %>%
            .[rsID == IndSigSNP] %>%
            .[nearestGene %in% genes[,symbol]] %>%
          .[, exonic := 0] %>%
          .[, CADD20 := 0] %>%
          setcolorder(., names(inGene))
          
        return(sentGene)    
    }

}) %>%
rbindlist %>%
setnames(., "nearestGene" , "gene") %>%
setnames(., "IndSigSNP", "sentinel") %>%
setnames(., "rsID", "mappedSNP") %>%
.[, uniqID := NULL]

fwrite(sentsAnnot, file = paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/",analysis,"_positionalMapping.csv"), sep = ",")



eQTLsWide <- allEQTLs  %>%
   .[ , c("gene", "sentinel", "eQTLtissue")] %>%
   .[, tissue := case_when(eQTLtissue == "EyeGEx" ~ "eQTLRetina",
            eQTLtissue %in% c("Cells_EBV-transformed_lymphocytes", "Whole_Blood") ~ "eQTLBlood",
            eQTLtissue %like% "Brain%" ~ "eQTLBrain" )] %>%
    dcast(., gene + sentinel  ~ tissue, fun.aggregate = function(x) as.numeric(length(x) > 0) )

ciWide <-  ci %>%
   .[ , c("gene", "sentinel", "tissue")] %>%
    .[, tissue := "chromatinInteractionCortex"] %>%
    dcast(., gene + sentinel  ~ tissue, fun.aggregate = function(x) as.numeric(length(x) > 0) )


allAnnot <- sentsAnnot %>%
    .[, .(sentinel, gene, exonic, CADD20)] %>%
    .[, proximity := 1] %>%
    unique %>%
    full_join(., eQTLsWide, by = c("sentinel", "gene")) %>%
    full_join(.,ciWide, by = c("sentinel", "gene"))
 allAnnot[is.na(allAnnot)] <- 0 

cols <-  c("sentinel", "gene", "proximity", "exonic", "CADD20", "eQTLBlood", 
           "eQTLBrain", "eQTLRetina","chromatinInteractionCortex")
missingCols <- setdiff(cols, names(allAnnot))

print(paste(analysis,"missing columns:", missingCols))
# add missing column with 0 values
if(length(missingCols) > 0) {
  allAnnot <- allAnnot[, (missingCols) := 0] %>%
    setcolorder(., cols)
}

allAnnot <- allAnnot[, sum_cols := rowSums(.SD), .SDcols = c("proximity", "exonic", "CADD20", "eQTLBlood", "eQTLBrain", "eQTLRetina", "chromatinInteractionCortex")] %>%
  .[, analysis := paste0(analysis) ]

return(allAnnot)

fwrite(allAnnot, file = paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/",analysis,"_summaryAnnotations.csv"), sep = ",")

}) %>%
  rbindlist



fpcSNPs <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/FPCall/FUMA_job248242/snps.txt") %>%
  .[rsID==IndSigSNP] %>%
  .[, sentinel := rsID] %>%
  .[, .(sentinel, chr, pos)]
pixelSNPs <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/allPixels/FUMA_job248074/snps.txt") %>%
  .[rsID==IndSigSNP] %>%
  .[, sentinel := rsID] %>%
  .[, .(sentinel, chr, pos)]

allSNPs <- rbind(fpcSNPs, pixelSNPs) %>%
  unique

# annotPOS <- fpcSNPs[annot, on = "sentinel"] %>%
#               .[pixelSNPs, on = c("sentinel", "chr", "pos")]


annotPOS <- allSNPs[annot, on = "sentinel"] %>%
  .[analysis := NULL] 


fpcResults <- lapply(c(1:22, "X"), function(x) {
  chrNo <- ifelse(x=="X", 23, x)
  
  print(paste("chr",x))
  
  chrPOS <- annotPOS[chr==chrNo]
  
  if(nrow(chrPOS) > 0) {
    
  fpcRes <- lapply(c(1:6), function(i) {
    print(paste(i))
    
    paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/results/chr",x,"/chr",x,"EUR.fpc",i,".glm.linear") %>%
      fread(., select = c("ID", "POS", "BETA", "P")) %>%
      setnames(., c("sentinel", "pos", "beta", "P")) %>%
      .[, chr := chrNo]
      .[, FPC := paste0("FPC",i)] %>%
      .[pos %in% chrPOS[,pos]] %>%
      return(.)
  }) %>%
    rbindlist %>%
    .[, Padj := P / 6] %>%
    dcast(., ...  ~ FPC, value.var = c("beta", "P", "Padj"))
  
  } else{
    print(paste("no SNPs for chr",x))
    return(NULL)
  }
  
}) %>%
  rbindlist



# fpcResults <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/sentinels/allChr_sentinel_clumpThresh0.001_withOverlap.txt", select = c("chr", "POS", "ID", "REF", "ALT", "FPC", "P", "clumpFPCs")) %>%
#   setnames(., "POS", "pos") %>%
#   fpcSNPs[., on = c("chr", "pos")]



pixelResults <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/allSentinelsAllPixelsResults_clumpThresh0.001_withOverlap.csv") %>%
  .[, Padj := P/29041] %>%
  .[, .SD[which.min(Padj)], by = ID] %>%
  .[, .(chrom, POS, BETA, P, Padj)] %>%
  setnames(., c("chr", "pos", "beta_allPixels", "P_allPixels", "Padj_allPixels"))


annotResults <- fpcResults %>%
 .[,analysis := NULL] %>%
  unique %>%
  pixelResults[., on = c("chr", "pos")] %>%
  .[, .SD[which.max(sum_cols)], by = sentinel]

analysesCols <- names(annotResults)[names(annotResults) %like any% c("Padj_%", "beta_%")]
betaCols <- names(annotResults)[names(annotResults) %like% "beta_%"]
PadjCols <- names(annotResults)[names(annotResults) %like% "Padj_%"]

betasLong <- annotResults %>%
  reshape2::melt(., id.vars = c("sentinel", "gene"), 
                 measure.vars = betaCols,
                 variable.name = "analysis", value.name = "beta") %>%
  as.data.table %>%
  .[, analysis := str_split(analysis, "_", simplify = TRUE)[, 2]]

pAdjLong <- annotResults %>%
  reshape2::melt(., id.vars = c("sentinel", "gene"), 
                 measure.vars = PadjCols,
                 variable.name = "analysis", value.name = "Padj") %>%
  as.data.table %>%
  .[, analysis := str_split(analysis, "_", simplify = TRUE)[, 2]]

plotAnalyses <-  betasLong[pAdjLong, on = c("sentinel", "gene", "analysis")] %>%
  .[, plotVal := sign(beta)*log(Padj, 10)]  %>%
  .[, .(gene, analysis, beta, Padj, plotVal)]

geneOrder <- plotAnalyses %>%
  .[, .SD[which.min(Padj)], by = gene] %>%
  .[order(Padj)] %>%
  .[,gene]



colours <- rev(brewer.pal(9,"RdBu"))
lim <- plotAnalyses[,plotVal] %>% abs %>% max(., na.rm=T) %>% ceiling()
  
ggplot(plotAnalyses, aes(x = analysis, y = gene, fill = plotVal)) +
  geom_tile() +
#  scale_fill_manual(values = color_scale) +
  labs(y = "Genes", title = "Implicated Genes") +
#  facet_wrap(. ~ cat, scales = "free_x" ) +
#  theme(legend.position = "none") +
  scale_fill_gradientn(colors=colours, limits = c(-lim, lim), breaks = c(-lim, -40, -30, -20, -10, 0, 10, 20, 30, 40, lim) ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_discrete(limits = rev(geneOrder) )





# convert data frame to long format using melt
plotAnnot1 <- melt(annotResults, id.vars = c("gene"), measure.vars = c("proximity", "eQTLBlood", "eQTLBrain", "eQTLRetina", "chromatinInteractionCortex")) %>%
  .[, .(sum_cols = sum(sum_cols), value = sum(value)), by = .(gene, variable)] %>%
  .[, plotVal := ifelse(value > 0, 1, 0)] %>%
  .[order(sum_cols, decreasing = T)] %>%
  .[, cat := "annotation"] %>%
  .[, .(gene, variable, sum_cols, cat, plotVal)]

plotAnnot2 <- dcast(annot, gene ~ analysis, fun.aggregate = function(x) as.numeric(length(x) > 0) ) %>%
  .[, sum_cols := rowSums(.SD), .SDcols = c("pixelWise", "FPC1", "FPC2", "FPC3", "FPC4", "FPC5", "FPC")] %>%
  melt(., id.vars = c("gene", "sum_cols")) %>%
  .[, plotVal := ifelse(value > 0, 1, 0)] %>%
  .[order(sum_cols, decreasing = T)] %>%
  .[, cat := "analysis"] %>%
  .[, .(gene, variable, sum_cols, cat, plotVal)]

orderList <- c("FPC1", "FPC2", "FPC3", "FPC4", "FPC5", "FPC6", "pixelWise", "proximity", "eQTLRetina", "eQTLBlood", "eQTLBrain", "chromatinInteractionCortex")

plotAnnot <- rbind(plotAnnot1, plotAnnot2) %>%
  arrange(factor(variable, levels = orderList), desc(plotVal))

keepGene <- plotAnnot[sum_cols>3,gene] %>% unique


plotAnnotA <- plotAnnot[gene %in% keepGene[1: (length(keepGene)/2)]] %>%
  arrange(factor(variable, levels = orderList), desc(plotVal))
plotAnnotB <- plotAnnot[gene %in% keepGene[-(1: (length(keepGene)/2))]] %>%
  arrange(factor(variable, levels = orderList), desc(plotVal))

geneOrderA <- unique(plotAnnotA[,gene])  %>%
  rev
geneOrderB <-  unique(plotAnnotB[,gene])  %>%
  rev

# create binary color scale
color_scale <- c("white", "blue")

# create heat map using geom_tile
ggplot(plotAnnotA, aes(x = variable, y = gene, fill = as.factor(plotVal))) +
  geom_tile() +
  scale_fill_manual(values = color_scale) +
  labs(y = "Genes", title = "Implicated Genes") +
  facet_wrap(. ~ cat, scales = "free_x" ) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_discrete(limits = geneOrderA)

ggplot(plotAnnotB, aes(x = variable, y = gene, fill = as.factor(plotVal))) +
  geom_tile() +
  scale_fill_manual(values = color_scale) +
  labs(y = "Genes", title = "Implicated Genes") +
  facet_wrap(. ~ cat, scales = "free_x" ) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_discrete(limits = geneOrderB)























# convert data frame to long format using melt
plotAnnot1 <- melt(annot, id.vars = c("gene","sum_cols"), measure.vars = c("proximity", "eQTLBlood", "eQTLBrain", "eQTLRetina", "chromatinInteractionCortex")) %>%
  .[, .(sum_cols = sum(sum_cols), value = sum(value)), by = .(gene, variable)] %>%
  .[, plotVal := ifelse(value > 0, 1, 0)] %>%
  .[order(sum_cols, decreasing = T)] %>%
  .[, cat := "annotation"] %>%
  .[, .(gene, variable, sum_cols, cat, plotVal)]

plotAnnot2 <- dcast(annot, gene ~ analysis, fun.aggregate = function(x) as.numeric(length(x) > 0) ) %>%
  .[, sum_cols := rowSums(.SD), .SDcols = c("pixelWise", "FPC1", "FPC2", "FPC3", "FPC4", "FPC5", "FPC")] %>%
  melt(., id.vars = c("gene", "sum_cols")) %>%
  .[, plotVal := ifelse(value > 0, 1, 0)] %>%
  .[order(sum_cols, decreasing = T)] %>%
  .[, cat := "analysis"] %>%
  .[, .(gene, variable, sum_cols, cat, plotVal)]

orderList <- c("FPC1", "FPC2", "FPC3", "FPC4", "FPC5", "FPC6", "pixelWise", "proximity", "eQTLRetina", "eQTLBlood", "eQTLBrain", "chromatinInteractionCortex")

plotAnnot <- rbind(plotAnnot1, plotAnnot2) %>%
  arrange(factor(variable, levels = orderList), desc(plotVal))
  
keepGene <- plotAnnot[sum_cols>3,gene] %>% unique
  




plotAnnotA <- plotAnnot[gene %in% keepGene[1: (length(keepGene)/2)]] %>%
  arrange(factor(variable, levels = orderList), desc(plotVal))
plotAnnotB <- plotAnnot[gene %in% keepGene[-(1: (length(keepGene)/2))]] %>%
  arrange(factor(variable, levels = orderList), desc(plotVal))

geneOrderA <- unique(plotAnnotA[,gene])  %>%
  rev
geneOrderB <-  unique(plotAnnotB[,gene])  %>%
  rev

# create binary color scale
color_scale <- c("white", "blue")

# create heat map using geom_tile
ggplot(plotAnnotA, aes(x = variable, y = gene, fill = as.factor(plotVal))) +
  geom_tile() +
  scale_fill_manual(values = color_scale) +
  labs(y = "Genes", title = "Implicated Genes") +
  facet_wrap(. ~ cat, scales = "free_x" ) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_discrete(limits = geneOrderA)

ggplot(plotAnnotB, aes(x = variable, y = gene, fill = as.factor(plotVal))) +
  geom_tile() +
  scale_fill_manual(values = color_scale) +
  labs(y = "Genes", title = "Implicated Genes") +
  facet_wrap(. ~ cat, scales = "free_x" ) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_discrete(limits = geneOrderB)





gwasCat <- fread(gwascatFile) %>%
  .[snps[, .(IndSigSNP, rsID, r2)] , on = c("IndSigSNP" = "IndSigSNP", "snp" = "rsID")] %>%
  .[!is.na(Trait)] %>%
  .[r2 > 0.5] %>%
  .[, .(IndSigSNP, snp, r2, PMID, Trait)] %>% 
  .[, .(PMIDs = paste(PMID, collapse = ";")), by = .(IndSigSNP, Trait)]


gwasPlot <- gwasCat[, .N, by = Trait] %>%
  .[N >1] %>%
  .[order(N, decreasing = T)] %>%
  .[, Trait := fct_reorder(Trait, N,  .desc = FALSE)]

# create bar plot of trait counts

ggplot(gwasPlot, aes(y = Trait, x = N)) +
  geom_bar(stat = "identity") +
  labs(y = "Trait", x = "Count", title = "GWAS catalogue Traits")

## FUMA files
#  params.config
#  IndSigSNPs.txt
#  leadSNPs.txt
#  GenomicRiskLoci.txt
#  snps.txt
#  ld.txt
#  annov.txt
#  annov.stats.txt
#  annot.txt
#  genes.txt
#  eqtl.txt
#  ci.txt
#  ciSNPs.txt
#  ciProm.txt
#  gwascatalog.txt
