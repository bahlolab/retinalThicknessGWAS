#!/usr/bin/env Rscript


library(data.table)
library(magrittr)
library(tidyverse)
library(optparse)
library(DescTools)
library(stringr)
library(RColorBrewer)
library(ggside)
require(gridExtra)
library(patchwork)

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


pixelWiseSentinels <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/rsIDsSentinelsPixelwiseBonferroniSig.txt")

fpcSentinels <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/rsIDsSentinelsFPCsBonferroniSig.txt")

mousepheno <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/abnormal_eye_morphology.tsv")
omim <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/OMIM-Gene-Map-Retrieval.tsv", skip = 4, sep = "\t", fill=T)

mouseGenes <- mousepheno[,Gene] %>% toupper() %>% unique
omimGenes <-  omim[, `Approved Symbol`] %>% unique 

# annot <- lapply(c("pixelWise", "FPC1", "FPC2", "FPC3", "FPC4", "FPC5", "FPC6"), function(analysis) {
annot <- lapply(c("pixelWise", "FPCall"), function(analysis) {
    
  if(analysis=="pixelWise") { 
    dir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/allPixels/FUMA_job248074"
    keep <- pixelWiseSentinels[,rsID]}

  if(analysis=="FPCall") { 
    dir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/FPCall/FUMA_job248242"
    keep <- fpcSentinels[,rsID]}
  
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
snps <- fread(snpsFile) %>%
  .[IndSigSNP %in% keep]

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

allAnnot <- allAnnot %>%
  .[, OMIM := ifelse(gene %in% omimGenes, 1, 0)] %>%
  .[, mousePhenotype := ifelse(gene %in% mouseGenes, 1, 0)] %>%
  .[, sum_cols := rowSums(.SD), .SDcols = c("proximity", "exonic", "CADD20", "eQTLBlood", "eQTLBrain", "eQTLRetina", "chromatinInteractionCortex", "OMIM", "mousePhenotype")] %>%
  .[, analysis := paste0(analysis) ]


fwrite(allAnnot, file = paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/",analysis,"_summaryAnnotations.csv"), sep = ",")

return(allAnnot)

}) %>%
  rbindlist 

setcolorder(annot, c("analysis", "sentinel", "gene", "proximity", "exonic",  "CADD20",
                     "eQTLRetina", "eQTLBlood",  "eQTLBrain", "chromatinInteractionCortex",
                     "OMIM", "mousePhenotype"))
fwrite(annot, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/geneSummaryAnnotations.csv", sep = ",")


## GWAS catalogue annotations

lapply(c("pixelWise","FPCall", "FPC1", "FPC2", "FPC3", "FPC4", "FPC5", "FPC6"), function(analysis) {

    if(analysis=="pixelWise") { 
      dir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/allPixels/FUMA_job248074"
      keep <- pixelWiseSentinels[,rsID]}
    
    if(analysis=="FPCall") { 
      dir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/FPCall/FUMA_job248242"
      keep <- fpcSentinels[,rsID]}
    
    if(analysis=="FPC1") { 
      dir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/FPC1/FUMA_job246508"
      keep <- fpcSentinels[,rsID]}

    if(analysis=="FPC2") { 
      dir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/FPC2/FUMA_job246510"
      keep <- fpcSentinels[,rsID]}

    if(analysis=="FPC3") {
      dir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/FPC3"
      keep <- fpcSentinels[,rsID]}

    if(analysis=="FPC4") {
      dir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/FPC4/FUMA_job246513"
      keep <- fpcSentinels[,rsID]}

    if(analysis=="FPC5") {
      dir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/FPC5/FUMA_job246514"
      keep <- fpcSentinels[,rsID]}

    if(analysis=="FPC6") {
      dir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/FPC6/FUMA_job246516"
      keep <- fpcSentinels[,rsID]}

    gwascatFile <- paste0(dir,"/gwascatalog.txt")
    snpsFile <- paste0(dir,"/snps.txt")
    
    snps <- fread(snpsFile) %>%
      .[IndSigSNP %in% keep]
    


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

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/",analysis,"_gwasCat.png"), width = 600, height = 600)
print({ ggplot(gwasPlot[1:20], aes(y = Trait, x = N)) +
  geom_bar(stat = "identity") +
  labs(y = "Trait", x = "Count", title = "GWAS catalogue Traits") 
})
dev.off()

})
  
  
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
