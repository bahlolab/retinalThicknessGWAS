library(data.table)
library(magrittr)
library(tidyverse)
library(topr)
library(ieugwasr)


## function - modified version of topr package. uses topr data
annotateGene <- function(variants, refGenes){
  if("POS" %in% colnames(variants) & "CHROM" %in% colnames(variants)){
    if(length(variants$POS) > 1000){
      print(paste("The dataset includes [",length(variants$POS),"] variants. This may take a while...", sep=""))
    }
    for(i in seq_along(variants$POS)){
      if(length(variants$POS) > 1000){
        if(i %% 1000==0) {
          # Print on the screen some message
          print(paste(i," variants annotated", sep=""))
        }
      }

      nearest_gene <- NULL
      variant <- variants[i,]
      chr <- gsub("chr", "", variant$CHROM)
      chr <- paste("chr",chr,sep="")
    
      genes_on_chr <- refGenes
      
 
      within_gene <-  genes_on_chr %>% dplyr::filter(gene_end >= variant$POS & gene_start <= variant$POS)
     
      if(length(within_gene$gene_symbol) > 0 ){  #TODO: order the genes by their biotype, and pull out the top one
        if(length(within_gene) == 1){ nearest_gene <- within_gene$gene_symbol }
        else{
          prot_coding <- within_gene %>% dplyr::filter(biotype=="protein_coding")
          if(length(prot_coding$gene_symbol) > 0){  nearest_gene <- prot_coding %>% utils::head(n=1)}
          else{  nearest_gene <- within_gene %>% utils::head(n=1) }
        }
      }

      if(! is.null(nearest_gene)){
        variants[i,"Gene_Symbol"] <- nearest_gene$gene_symbol
      }else{
        variants[i,"Gene_Symbol"] <- "intergenic"
      }
    }
  }
  else{
    stop("Cannot find the columns CHROM and POS in the input data. Add the required columns and try again, or rename existing columns, e.g. df=df %>% dplyr::rename(CHROM=yourColname)")
  }
  return(variants)
}


rsids <- fread("/vast/scratch/users/jackson.v/retThickness/GWAS/tmp/chr22.vcf", header=F) %>%
setnames(., c("CHROM", "POS", "rsid", "REF", "ALT"))


variants <- fread("/vast/scratch/users/jackson.v/retThickness/GWAS/annot/chr22.afreq") %>%
setnames(., "#CHROM", "CHROM")

freqs <- fread("/vast/scratch/users/jackson.v/retThickness/GWAS/annot/chr22.afreq") %>%
setnames(., "#CHROM", "CHROM")

variants <- fread("/vast/scratch/users/jackson.v/retThickness/GWAS/annot/chr22.pvar") %>%
setnames(., "#CHROM", "CHROM") %>%
freqs[., on= c("CHROM", "ID", "REF", "ALT")] %>%
rsids[., on = c("CHROM", "POS", "REF", "ALT")] %>%
.[, SNPID := case_when(ID %like% "rs" ~ ID, 
                        is.na(rsid) ~ ID,
                        T ~ rsid)] %>%
.[, .(SNPID, CHROM, POS, REF, ALT, ALT_FREQS)] %>%
unique


refGenes <- toprdata::ENSGENES_37 %>% dplyr::filter(chrom == 22) %>% dplyr::arrange(gene_start)
variantsAnnot <- annotateGene(variants, refGenes)

fwrite(variantsAnnot, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/finalResultsEUR/chr22SNPinfo.txt")
