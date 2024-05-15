#!/usr/bin/env Rscript

library(optparse)
# BiocManager::install(c("seqarray", "VariantAnnotation", "SeqVarTools"))
# Load required packages
library(SeqArray)
library(VariantAnnotation)
library(SeqVarTools)
library(VariantAnnotation)
library(data.table)
library(magrittr)
library(dplyr)

# Parse command line arguments
option_list <-  list(
  make_option(c("-c", "--chr"), type="integer", default=NULL,
              help="chromosome number", metavar="integer"));

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

chr <-  opt$chr

setwd("/vast/scratch/users/jackson.v/harmonisingData")

  
  print(paste("Processing chromosome", chr))
  
  # Read PLINK files
  plinkData <- fread(paste0("/vast/scratch/users/jackson.v/retThickness/GWAS/geneticData/cleanedEURData/plinkBin/EUR_minMaf0.005_minInfo0.8_chr",chr,".bim"), select = c(1:2,4:6)) %>%
    setnames(., c("CHR", "IDukbb", "POS", "A1ukbb", "A2ukbb"))
  
  # Read VCF file
  vcfFile <- paste0("/vast/scratch/users/jackson.v/harmonisingData/inVcf/chr",chr,".dose.vcf.gz")
  
  ## create gds
  seqVCF2GDS(vcfFile , paste0("./tmp/chr",chr,".gds"), storage.option="ZIP_RA")
  
  ## open gds
  vcfData <- seqOpen( paste0("./tmp/chr",chr,".gds"))
  
  ## extract SNP info
  vcfVars <- data.table(idxeqtl = seqGetData(vcfData, "variant.id"),
                        IDeqtl = seqGetData(vcfData,"annotation/id"), 
                        CHR = seqGetData(vcfData, "chromosome") %>% as.integer,
                        POS = seqGetData(vcfData, "position"),
                        alleles = seqGetData(vcfData, "allele")) %>%
    .[, c("A1eqtl", "A2eqtl") := tstrsplit(alleles, ",")]
  
  
  ## merge UK Biobank plink file and vcf
  merged <- plinkData[vcfVars, on = c("CHR", "POS")]
  
  ## subset to variants present in ukbb
  shared <- merged[!is.na(IDukbb)]
  
  ## check alleles
  print(paste("There are", shared %>% nrow, "shared variants in total"))
  print(paste("There are", shared[A1ukbb==A2eqtl & A2ukbb==A1eqtl] %>% nrow, "variants with matching alleles"))
  print(paste("There are", shared[!(A1ukbb==A2eqtl & A2ukbb==A1eqtl)] %>% nrow,"mismatches"))
  
  mismatches <- shared[!(A1ukbb==A2eqtl & A2ukbb==A1eqtl)] 
  
  ## Identify non A/T/C/G variants
  nonSNPs <- mismatches[!(A1ukbb %in% c("A", "C", "G", "T") & A2ukbb %in% c("A", "C", "G", "T"))]
  print(paste(nrow(nonSNPs), "non-SNPs to be removed"))
  
  snpsToResolve <-  mismatches[A1ukbb %in% c("A", "C", "G", "T") & A2ukbb %in% c("A", "C", "G", "T")] %>%
                    .[, A1flip := case_when(A1ukbb=="A" & A2eqtl=="T" ~ "flip",
                                            A1ukbb=="C" & A2eqtl=="G" ~ "flip",
                                            A1ukbb=="G" & A2eqtl=="C" ~ "flip",
                                            A1ukbb=="T" & A2eqtl=="A" ~ "flip",
                                            T ~ "noMatch")]  %>%
                    .[, A2flip := case_when(A2ukbb=="A" & A1eqtl=="T" ~ "flip",
                                            A2ukbb=="C" & A1eqtl=="G" ~ "flip",
                                            A2ukbb=="G" & A1eqtl=="C" ~ "flip",
                                            A2ukbb=="T" & A1eqtl=="A" ~ "flip",
                                            T ~ "noMatch")] %>%
                    .[, snpFlip := case_when(A1flip == "flip" & A2flip == "flip" ~ "flip",
                                             T ~ "noMatch")]
  
  ## check remaining are true allele mismatches
  print(paste(nrow(snpsToResolve[snpFlip == "noMatch"]), "out of", nrow(snpsToResolve), "variants have mismatched alleles and need to be removed"))
  print(paste(nrow(snpsToResolve[snpFlip == "flip"]), "variants could potentially be recovered as strand flips"))
  
  if(nrow(snpsToResolve[snpFlip == "flip"]) == 0) {
    
    ## obtain SNPs to be kept
    keep <- shared[A1ukbb==A2eqtl & A2ukbb==A1eqtl] 
    
    ## filter to SNPs to keep and output vcf
    seqSetFilter(vcfData, variant.id = keep[,idxeqtl])
    seqGDS2VCF(vcfData, paste0("./tmp/chr",chr,"_filtered.vcf"))
    seqClose(vcfData)
  
    # read vcf with VariantAnnotation package
    vcfFilt <- readVcf(paste0("./tmp/chr",chr,"_filtered.vcf"))
    
    ## double check IDs match
    print("Do all of the IDs match?")
    ifelse(identical(vcfFilt@rowRanges@ranges@NAMES, keep[,IDeqtl]), print("Yes!"), print("No!"))
    
    # Modify the variant IDs
    vcfFilt@rowRanges@ranges@NAMES <- keep[,IDukbb]
    
    # Write the modified VCF file
    writeVcf(vcfFilt, paste0("outVcf/chr",chr,"_harmonised.vcf"))           
  
    rm(list=ls())
    gc()
    
  }
