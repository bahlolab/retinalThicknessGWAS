library(magrittr)
library(data.table)
library(tidyverse)


sentinelUpload <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/pixelWiseFUMA.txt")
eqtl <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/eqtl.txt")
genes <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/genes.txt")
snps <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/snps.txt")

retinaEQTLs <- snps[, .(IndSigSNP, uniqID, r2)] %>%
    .[eqtl[db == "EyeGEx", .(uniqID, symbol, p, FDR)], on = "uniqID"] %>%
    .[, .SD[which.min(p)], by = c("symbol", "IndSigSNP")]

allEQTLs <- snps[, .(IndSigSNP, uniqID, r2)] %>%
    .[eqtl[, .(uniqID, tissue, symbol, p, FDR)], on = "uniqID"] %>%
    .[, .SD[which.min(p)], by = c("symbol", "IndSigSNP", "tissue")] %>%
    .[r2 > 0.5] %>%
    setnames(., c("gene", "sentinel", "eQTLtissue", "eSNP", "eSNPr2", "P", "FDR")) %>%
    .[gene %in% genes[,symbol]]


ciSNPs <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/ciSNPs.txt") %>%
.[snps[, .(IndSigSNP, uniqID, r2)] , on = "uniqID"] %>%
.[!is.na(rsID)] %>%
.[, .(rsID, IndSigSNP, uniqID, r2)] %>%
unique

ci <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/ci.txt") %>%
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
    setnames("IndSigSNP", "sentinel") %>%
    .[gene %in% genes[,symbol]]


position <- snps[r2 > 0.5] %>%
    .[, .(uniqID, rsID, r2, IndSigSNP, nearestGene, dist, func)]


sents <- position[,IndSigSNP] %>% unique

sentsAnnot <- lapply(sents, function(s) {

    print(paste(s))

    sentsMaps <- position[IndSigSNP == s]

    inGene <- sentsMaps %>%
        tidyr::separate_rows(., nearestGene, dist, sep = ":") %>%
        as.data.table() %>%
        .[dist==0] %>%
        .[, .SD[which.max(r2)], by = c("IndSigSNP", "nearestGene")] %>%
        .[nearestGene %in% genes[,symbol]]


    if(nrow(inGene) > 0) {
        return(inGene)
    } else {
        sentGene <- sentsMaps %>%
            tidyr::separate_rows(., nearestGene, dist, sep = ":") %>%
            as.data.table() %>%
            .[rsID == IndSigSNP] %>%
            .[nearestGene %in% genes[,symbol]] %>%
          setcolorder(., names(inGene))
          
        return(sentGene)    
    }

}) %>%
rbindlist %>%
setnames(., "nearestGene" , "gene") %>%
setnames(., "IndSigSNP", "sentinel") %>%
setnames(., "rsID", "mappedSNP") %>%
.[, uniqID := NULL]




eQTLsWide <- allEQTLs  %>%
   .[ , c("gene", "sentinel", "eQTLtissue")] %>%
   .[, tissue := case_when(eQTLtissue == "EyeGEx" ~ "eQTLRetina",
            eQTLtissue %in% c("Cells_EBV-transformed_lymphocytes", "Whole_Blood") ~ "eQTLBlood",
            eQTLtissue %like% "Brain" ~ "eQTLBrain" )] %>%
    dcast(., gene + sentinel  ~ tissue, fun.aggregate = function(x) as.numeric(length(x) > 0) )

ciWide <-  ci %>%
   .[ , c("gene", "sentinel", "tissue")] %>%
    .[, tissue := "chromatinInteractionCortex"] %>%
    dcast(., gene + sentinel  ~ tissue, fun.aggregate = function(x) as.numeric(length(x) > 0) )


allAnnot <- sentsAnnot %>%
    .[, .(sentinel, gene)] %>%
    .[, positionalMapping := 1] %>%
    unique %>%
    full_join(., eQTLsWide, by = c("sentinel", "gene")) %>%
    full_join(.,ciWide, by = c("sentinel", "gene"))
 allAnnot[is.na(allAnnot)] <- 0 

allAnnot <- allAnnot[, sum_cols := rowSums(.SD), .SDcols = c("positionalMapping", "eQTLBlood", "eQTLBrain", "eQTLRetina", "chromatinInteractionCortex")]

gwasCat <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/data/gwascatalog.txt") %>%
  .[snps[, .(IndSigSNP, rsID, r2)] , on = c("IndSigSNP" = "IndSigSNP", "snp" = "rsID")] %>%
  .[!is.na(Trait)] %>%
  .[r2 > 0.5] %>%
  .[, .(IndSigSNP, snp, r2, PMID, Trait)] %>% 
  .[, .(PMIDs = paste(PMID, collapse = ";")), by = .(IndSigSNP, Trait)]



# convert data frame to long format using melt
plotAnnot <- melt(allAnnot, id.vars = c("sentinel", "gene", "sum_cols")) %>%
  .[order(sum_cols, decreasing = T)]

# create binary color scale
color_scale <- c("white", "blue")

# create heat map using geom_tile
ggplot(plotAnnot[sum_cols > 2], aes(x = variable, y = gene, fill = as.factor(value))) +
  geom_tile() +
  scale_fill_manual(values = color_scale) +
  labs(x = "Evidence", y = "Genes", title = "Implicated Genes")


gwasPlot <- gwasCat[, .N, by = Trait] %>%
  .[N > 5] %>%
  .[order(N, decreasing = T)] %>%
  .[, Trait := fct_reorder(Trait, N,  .desc = FALSE)]

# create bar plot of trait counts

ggplot(gwasPlot, aes(y = Trait, x = N)) +
  geom_bar(stat = "identity") +
  labs(y = "Trait", x = "Count", title = "GWAS catalogue Traits with more than five occurrence")

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
