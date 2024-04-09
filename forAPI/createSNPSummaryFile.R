#!/usr/bin/env Rscript

library(data.table)
library(magrittr)
library(tidyverse)
library(topr)
library(optparse)
# library(ieugwasr)

option_list <-  list(
 make_option(c("-c", "--chr"), type="integer", default=NULL,
             help="chromosome", metavar="integer"));

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
 
chromosome <-  ifelse(opt$chr == 23, "X", opt$chr)


## function - modified version of topr package. uses topr data
annotateGene <- function(variants, refGenes, chromosome){
  if("POS" %in% colnames(variants) & "CHROM" %in% colnames(variants)){
    if(length(variants$POS) > 1000){
      print(paste("The dataset includes [",length(variants$POS),"] variants. This may take a while...", sep=""))
    }
    
    variantsOut <- lapply(c(1:length(variants$POS)), function(i) {
      
      
      if(length(variants$POS) > 1000){
        if(i %% 1000==0) {
          # Print on the screen some message
          print(paste(i," variants annotated", sep=""))
        }
      }
      
      nearest_gene <- NULL
      variant <- variants[i,]
      
      chr <- paste0("chr",chromosome)
      
      genes_on_chr <- refGenes
      
      within_gene <-  genes_on_chr %>% dplyr::filter(gene_end >= variant$POS & gene_start <= variant$POS)
      
      if(length(within_gene$gene_symbol) > 0 ){  #TODO: order the genes by their biotype, and pull out the top one
        if(length(within_gene) == 1){ nearest_gene <- within_gene }
        else{
          prot_coding <- within_gene %>% dplyr::filter(biotype=="protein_coding")
          if(length(prot_coding$gene_symbol) > 0){  nearest_gene <- prot_coding %>% utils::head(n=1)}
          else{  nearest_gene <- within_gene %>% utils::head(n=1) }
        }
      }
      
      if(! is.null(nearest_gene)){
        variant  <- variant[,Gene_Symbol := nearest_gene$gene_symbol]
      } else{
        variant  <- variant[,Gene_Symbol := "intergenic"]
      }
      return(variant)
    }) %>%
      rbindlist
    return(variantsOut)
  }
  else{
    stop("Cannot find the columns CHROM and POS in the input data. Add the required columns and try again, or rename existing columns, e.g. df=df %>% dplyr::rename(CHROM=yourColname)")
  }

}


# rsids <- paste0("/vast/scratch/users/jackson.v/retThickness/GWAS/annot/chr",chr,".vcf") %>%
#   fread(., header=F) %>%
# setnames(., c("CHROM", "POS", "rsid", "REF", "ALT"))


result <-fread(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/results/chr",chromosome,"/chr",chromosome,"EUR.fpc1.glm.linear")) %>%
  setnames(., "#CHROM", "CHROM") %>%
  .[, .(CHROM, POS, ID, REF, ALT, A1, A1_FREQ)]

variants <- fread(paste0("/vast/scratch/users/jackson.v/retThickness/GWAS/annot/chr",chromosome,".pvar")) %>%
setnames(., "#CHROM", "CHROM") %>%
result[., on= c("CHROM", "ID", "REF", "ALT")] %>%
# rsids[., on = c("CHROM", "POS", "REF", "ALT")] %>%
# .[, SNPID := case_when(ID %like% "rs" ~ ID, 
#                         is.na(rsid) ~ ID,
#                         T ~ rsid)] %>%
  .[, SNPID := ID] %>%
  .[, effAllele := A1] %>%
  .[, noneffAllele := case_when(A1 == REF ~ ALT,
                                A1 == ALT ~ REF)] %>%
  .[, effAlleleFreq := A1_FREQ] %>%
  .[, .(SNPID, CHROM, POS, effAllele, noneffAllele, effAlleleFreq)] %>%
unique


refGenes <- toprdata::ENSGENES_37 %>% dplyr::filter(chrom == paste0("chr",chromosome)) %>% dplyr::arrange(gene_start)
variantsAnnot <- annotateGene(variants, refGenes, chromosome)

fwrite(variantsAnnot, file = paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/finalResultsEUR/chr",chromosome,"SNPinfo.txt"))
