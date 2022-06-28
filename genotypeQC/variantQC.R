#!/usr/bin/env Rscript

library(data.table)
library(magrittr)
library(tidyverse)
library(here)
library(optparse)

option_list <-  list(
  make_option(c("-s", "--snpQC"), type="character", default=NULL,
              help="file containing SNP QC info (ukb_snp_qc.txt)", metavar="character"),
  make_option(c("-m", "--mfiDir"), type="character", default=NULL,
              help="directory where the per chromosome UKBiobank MAF/INFO files (ukb_mfi_chr*_v3.txt) are located", metavar="character"),
	make_option(c("-o", "--outDir"), type="character", default=NULL,
              help="output directory", metavar="character")
);

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

snpQcFile <- opt$snpQC
mfiDir <- opt$mfiDir
outDir <- opt$outDir

# snpQcFile <- "data/ukb_snp_qc.txt"
# mfiDir <- "data"
# outDir <- "qcFiles"
print(paste0("Reading in SNP MAF INFO data from ",mfiDir,"/ukb_mfi_chr*_v3.txt"))

# Read in MAFs and imputation quality
chr <- c(1:22,"X","XY")
infoFreq <- lapply(chr, function(x) {
      paste0("ukb_mfi_chr",x,"_v3.txt") %>%
      here(mfiDir,.) %>%
      fread(.) %>%
      setnames(., c("SNP", "rsid", "POS", "Allele1", "Allele2", "MAF", "MinorAllele", "INFO")) %>%
      .[, CHR := x] %>%
      .[, .(SNP, rsid, CHR, POS, MAF, INFO)]
    }) %>%
  rbindlist

print(paste("writing list of variants to include under different filters to",outDir))
## Different filters
filt2 <- infoFreq[MAF>=0.001 & INFO>=0.8, SNP]
filt3 <- infoFreq[MAF>=0.001 & INFO>=0.5, SNP]

write.table(filt2, file=here(outDir, "snpInclude_minMaf0.001_minInfo0.8.txt"), row.names=F, col.names=F, quote=F)
write.table(filt3, file=here(outDir, "snpInclude_minMaf0.001_minInfo0.5.txt"), row.names=F, col.names=F, quote=F)

## AltIDs
filt2rsid <- infoFreq[MAF>=0.001 & INFO>=0.8, rsid]
filt3rsid <- infoFreq[MAF>=0.001 & INFO>=0.5, rsid]

write.table(filt2rsid, file=here(outDir, "snpIncludeAltID_minMaf0.001_minInfo0.8.txt"), row.names=F, col.names=F, quote=F)
write.table(filt3rsid, file=here(outDir, "snpIncludeAltID_minMaf0.001_minInfo0.5.txt"), row.names=F, col.names=F, quote=F)

print(paste("writing list of variants to use for genetic relatedness matrix (grm) to",outDir))
# Identify suset of SNPs used in relatedness calculations - these to be used for glm.
# Use plink IDs here
snpQC <- fread(snpQcFile) %>%
  .[in_PCA==1, rs_id]

write.table(snpQC, file=here(outDir, "snpsInclude_grm_plink.txt"), row.names=F, col.names=F, quote=F)
