#!/usr/bin/env Rscript

library(data.table)
library(magrittr)
library(tidyverse)
library(here)
library(optparse)

option_list <-  list(
  make_option(c("-w", "--withdrawn"), type="character", default=NULL,
              help="withdrawn list", metavar="character"),
  make_option(c("-s", "--sampQC"), type="character", default=NULL,
              help="phenotype file, containing sample qc fields (with default UKBiobank codes as colum names)", metavar="character"),
	make_option(c("-o", "--outDir"), type="character", default=NULL,
              help="output directory", metavar="character")
);

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

withdrawnFile <-  opt$withdrawn
sampleQCFile <- opt$sampQC
outDir <- opt$outDir


# withdrawnFile <-  "data/withdrawn.csv"
# sampleQCFile <- "data/sampleQC.tab"
# outDir <- "qcFiles"

print(paste("Reading in Withdrawn IDs from",withdrawnFile))

withdrawnIDs <- fread(here(withdrawnFile)) %>%
  .[[1]]

print(paste("Reading in Sample QC info from",sampleQCFile))

sampleQC <- fread(here(sampleQCFile),
  select=c("f.eid", "f.31.0.0", "f.22001.0.0", "f.22028.0.0", "f.22029.0.0", "f.22030.0.0", "f.22006.0.0", "f.22027.0.0", "f.22021.0.0", "f.22019.0.0"), quote="")

setnames(sampleQC, c("eid", "sex", "geneticSex", "in_phasing_input_chr1_22", "in_phasing_input_chrx", "in_phasing_input_chrxy", "whiteBritish", "missingHetOutlier", "geneticKinship", "putative_sex_chromosome_aneuploidy"))


## read in related and unrelated lists

relatedWhiteBritIDs <- paste0(outDir,"/relatedWhiteBritIDs.txt") %>% fread
unrelatedWhiteBritIDs <- paste0(outDir,"/unrelatedWhiteBritIDs.txt") %>% fread
unrelatedNonWBIDs <- paste0(outDir,"/unrelatedNonWhiteBritIDs.txt") %>% fread


sampQC <- sampleQC %>%
  .[, withdrawn := case_when(eid %in% withdrawnIDs ~ 1,
                            T ~ 0)] %>%
  .[, sexMismatch := case_when(sex != geneticSex ~ 1,
                            T ~ 0)] %>%
  .[, excess_relatives := case_when(geneticKinship == 10 ~ 1,
                            T ~ 0)] %>%
  .[, exclude := case_when(in_phasing_input_chr1_22==0 ~ 1,
                           sexMismatch==1 ~ 1,
                           excess_relatives==1 ~ 1,
                           putative_sex_chromosome_aneuploidy==1 ~ 1,
                           withdrawn==1 ~ 1,
                           T ~ 0)]

sampExclusion <- sampQC %>%
  .[, .(eid, in_phasing_input_chr1_22, in_phasing_input_chrx, in_phasing_input_chrxy, whiteBritish, sexMismatch, excess_relatives, putative_sex_chromosome_aneuploidy, withdrawn, exclude)] %>%
  .[, relatedWhiteBrit := case_when(eid %in% relatedWhiteBritIDs[,ID] ~ 1,
    T ~ 0)] %>%
  .[, unrelatedWhiteBrit := case_when(eid %in% unrelatedWhiteBritIDs[,ID] ~ 1,
    T ~ 0)] %>%
  .[, unrelatedNonWBIDs  := case_when(eid %in% unrelatedNonWBIDs[,ID] ~ 1,
    T ~ 0)]


print(paste("Outputting lists to",outDir))

# write table of exclusions and why
write.table(sampExclusion, file=here(outDir, "/sampleGenoQC.csv"), row.names=F, col.names=T, quote=F, sep=",")

# write list of unrelated White British ids to include
paste("Writing list of",nrow(sampExclusion[exclude==0 & unrelatedWhiteBrit==1]), "unrelated, White British individuals")
cat("\n")
write.table(sampExclusion[exclude==0 & unrelatedWhiteBrit==1, eid], file=here(outDir, "/sampleIncludeUnrelatedWhiteBrit.txt"), row.names=F, col.names=F, quote=F)
write.table(sampExclusion[exclude==0 & unrelatedWhiteBrit==1, .(eid, eid)], file=here(outDir, "/sampleIncludeUnrelatedWhiteBrit_plink.txt"), row.names=F, col.names=F, quote=F)

# write list of related White British ids to include
paste("Writing list of",nrow(sampExclusion[exclude==0 & relatedWhiteBrit==1]), "related, White British individuals")
cat("\n")
write.table(sampExclusion[exclude==0 & relatedWhiteBrit==1, eid], file=here(outDir, "/sampleIncludeRelatedWhiteBrit.txt"), row.names=F, col.names=F, quote=F)
write.table(sampExclusion[exclude==0 & relatedWhiteBrit==1, .(eid, eid)], file=here(outDir, "/sampleIncludeRelatedWhiteBrit_plink.txt"), row.names=F, col.names=F, quote=F)

# write list of unrelated non-White British ids to include
paste("Writing list of",nrow(sampExclusion[exclude==0 &  unrelatedNonWBIDs==1]), "unrelated, non-White British individuals")
cat("\n")
write.table(sampExclusion[exclude==0 &  unrelatedNonWBIDs==1, eid], file=here(outDir, "/sampleIncludeUnrelatedNonWBIDs.txt"), row.names=F, col.names=F, quote=F)
write.table(sampExclusion[exclude==0 &  unrelatedNonWBIDs==1, .(eid, eid)], file=here(outDir, "/sampleIncludeUnrelatedNonWBIDs_plink.txt"), row.names=F, col.names=F, quote=F)
