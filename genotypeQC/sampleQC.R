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

dataDir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/genotypeQC/rawData/"

# withdrawnFile <-  "data/withdrawn.csv"
# sampleQCFile <- "data/sampleQC.tab"
# outDir <- "qcFiles"

print(paste("Reading in Withdrawn IDs from",withdrawnFile))

withdrawnIDs <- fread(here(withdrawnFile)) %>%
  .[[1]]

print(paste("Reading in Sample QC info from",sampleQCFile))

sampleQC <- fread(here(sampleQCFile),
  select=c("f.eid", "f.31.0.0", "f.22001.0.0", "f.22028.0.0", "f.22029.0.0", "f.22030.0.0",  "f.22027.0.0", "f.22021.0.0", "f.22019.0.0"), quote="")

setnames(sampleQC, c("eid", "sex", "geneticSex", "in_phasing_input_chr1_22", "in_phasing_input_chrx", "in_phasing_input_chrxy",  "missingHetOutlier", "geneticKinship", "putative_sex_chromosome_aneuploidy"))



sampGenoQC <- sampleQC %>%
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


excludeSamps <- sampGenoQC[exclude==1, eid]

ids <- paste0(outDir,"idLinkage.txt") %>% fread

sampQC  <- paste0(dataDir,"ukb32825.tab") %>%
  fread(., select = c(1, 1187, 1145:1154)) %>%
  setnames(. ,c("patIDhda", "related", paste0("PC",c(1:10))))

sampleFile <- paste0(dataDir,"cleanedEuro_chr1_chunk1.sample") %>% fread

relatedness <- paste0(dataDir,"ukb_rel_a36610_s488212.dat") %>% fread

ancestryLinkIDs <- paste0(dataDir,"ukb36610bridge31063.txt") %>% 
  fread %>%
  setnames(., c("patIDhda", "patIDlink"))

ancestry <- paste0(dataDir,"all_pops_non_eur_pruned_within_pop_pc_covs.tsv") %>% 
  fread %>%
  setnames(., "s", "patIDlink") %>%
  .[, .(patIDlink, pop)] %>%
  ancestryLinkIDs[., on = "patIDlink"]

## explore ancestry in full cohort (restricting to unrelated individuals)
sampQCPheno <- sampQC[ids, on = "patIDhda"] %>%
  .[ancestry, on = "patIDhda"] %>%
  .[ids, on = "patIDhda"] %>%
  .[!patIDhda %in% excludeSamps]


ggplot(sampQCPheno, aes(x = PC1, y = PC2, group = pop, color = pop)) +
  geom_point()

ggsave(file = paste0(outDir,"PC1vsPC2pop_withOCT.png"))

ggplot(sampQCPheno, aes(x = PC3, y = PC4, group = pop, color = pop)) +
  geom_point()

ggsave(file = paste0(outDir,"PC3vsPC4pop_withOCT.png"))


sampQCPheno[, pop] %>% table(., useNA = "always")
sampQCPheno[related==0, pop] %>% table(., useNA = "always")
sampQCPheno[related==1, pop] %>% table(., useNA = "always")

## look at ancestry of those not assigned to a super pop...
ggplot(sampQC[patIDhda %in% sampQCPheno[is.na(pop), patIDhda]], aes(x = PC1, y = PC2)) +
  geom_point()


## restrict to unrelated set
related <- relatedness[ID1 %in% sampQCPheno[!is.na(pop),patIDhda] & ID2 %in% sampQCPheno[!is.na(pop),patIDhda]]

relativesCount <- c(related[,ID1], related[,ID2]) %>%
  table %>%
  as.data.table %>%
  setnames(., c("ID", "N"))

multiRelatives <- relativesCount[N>1, ID]
multiRelatives %>% length

## number of individuals recovered to unrelated set
recovered <- related[(ID1 %in% multiRelatives | ID2 %in% multiRelatives)] %>%
  .[, .(ID1,ID2)] %>%
  unlist %>%
  .[!(. %in% multiRelatives)]

recovered %>% length

relatedFilt1 <- related[!(ID1 %in% multiRelatives | ID2 %in% multiRelatives)]

## double check all unique pairs
c(relatedFilt1[,ID1], relatedFilt1[,ID2]) %>% 
  table %>% 
  table

set.seed(4860)

keep <- rbinom(nrow(relatedFilt1), 1,.5)

removeRandom <- relatedFilt1[, remove := ifelse(keep==1, ID1, ID2)] %>%
  .[, remove]


finalSet <- sampQCPheno[!(patIDhda %in% multiRelatives | patIDhda %in% removeRandom)]

finalSet %>% .[,  .(pop)] %>% table(., useNA = "always")

excluded <- sampQCPheno[(patIDhda %in% multiRelatives | patIDhda %in% removeRandom)]

excluded %>% .[,  .(pop)] %>% table(., useNA = "always")

## check all accounted for
nrow(finalSet) + nrow(excluded)
nrow(sampQCPheno)

## write out sample lists
lapply(c("EUR", "CSA", "AFR"), function(anc) {
  
  ids <-  finalSet[pop==anc, patIDhda]
  
  print(paste(length(ids), "individuals with", anc, "ancestry."))
  idDT <- data.table(FID = ids,
                     IID = ids)
  
  fwrite(idDT, file = paste0(outDir,"sampleList_doubleIDs_",anc,".txt"), sep = "\t", na = "NA", quote = F)
  fwrite(idDT[,!"FID"], file = paste0(outDir,"sampleList_singleIDs_",anc,".txt"), sep = "\t", na = "NA", quote = F)
  
})


ancestryExclude <-  finalSet[!pop %in% c("EUR", "CSA", "AFR"), patIDhda]
EURset <-  finalSet[pop == "EUR", patIDhda]
CSAset <-  finalSet[pop == "CSA", patIDhda]
AFRset <-  finalSet[pop == "AFR", patIDhda]

sampExclusion <- sampGenoQC %>%
  .[, relatednessExclusion := case_when(eid %in% excluded ~ 1,
                                        T ~ 0)] %>%
  .[, ancestryExclusion := case_when(eid %in% ancestryExclude ~ 1,
                                        T ~ 0)]  %>%
  .[,finalEURset := case_when(eid %in% EURset ~ 1,
                                     T ~ 0)] %>%  
  .[,finalCSAset := case_when(eid %in% CSAset ~ 1,
                                    T ~ 0)] %>%
  .[,finalAFRset := case_when(eid %in% AFRset ~ 1,
                              T ~ 0)] 
  
# write table of exclusions and why
write.table(sampExclusion, file=here(outDir, "/sampleGenoQC.csv"), row.names=F, col.names=T, quote=F, sep=",")

