## Run using R/4.1.2

## Link mactel project and healthy dev aging project IDs -
## cleaned genetic data has latter ids

library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)


dataDir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/genotypeQC/rawData/"
outDir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/genotypeQC/output/"


ids <- paste0(outDir,"idLinkage.txt") %>% fread

sampQC  <- paste0(dataDir,"ukb32825.tab") %>%
  fread(., select = c(1, 1124, 1142, 1187, 1145:1154)) %>%
  setnames(. ,c("patIDhda", "selfRepEthnic", "whiteBrit", "related", paste0("PC",c(1:10))))

sampleFile <- paste0(dataDir,"cleanedEuro_chr1_chunk1.sample") %>% fread

relatedness <- paste0(dataDir,"ukb_rel_a36610_s488212.dat") %>% fread


## explore ancestry in full cohort (restricting to unrelated individuals)
sampQCPheno <- sampQC[ids, on = "patIDhda"] %>%
.[, selfRepEthnicGroup := case_when(selfRepEthnic %in% c(1, 1001:1003) ~ "White",
                                selfRepEthnic %in% c(2, 2001:2004) ~  "Mixed" ,
                              selfRepEthnic %in% c(3, 3001:3004) ~ "Asian/South Asian",
                            selfRepEthnic %in% c(4, 4001:4003) ~ "Black",
                          selfRepEthnic %in% 5 ~ "Chinese",
                        selfRepEthnic %in% 6 ~ "other")]

pcPlot <- ggplot(sampQCPheno, aes(x = PC1, y = PC2, group = whiteBrit, color = whiteBrit)) +
  geom_point()

ggsave(file = paste0(outDir,"PC1vsPC2whiteBrit.png"))

pcPlot <- ggplot(sampQCPheno, aes(x = PC1, y = PC2, group = selfRepEthnicGroup, color = selfRepEthnicGroup)) +
  geom_point()

ggsave(file = paste0(outDir,"PC1vsPC2selfReportEthnicity.png"))

sampQCPheno[, .(whiteBrit, selfRepEthnicGroup)] %>% table(., useNA = "always")
sampQCPheno[related==0, .(whiteBrit, selfRepEthnicGroup)] %>% table(., useNA = "always")
sampQCPheno[related==1, .(whiteBrit, selfRepEthnicGroup)] %>% table(., useNA = "always")


## restrict to unrelated set
related <- relatedness[ID1 %in% sampQCPheno[,patIDhda] & ID2 %in% sampQCPheno[,patIDhda]]

relativesCount <- c(related[,ID1], related[,ID2]) %>%
  table %>%
  as.data.table %>%
  setnames(., c("ID", "N"))

multiRelatives <- relativesCount[N>1, ID]
multiRelatives %>% length

## number of individuals recovered to unrelated set
recovered <- related[(ID1 %in% multiRelatives | ID2 %in% multiRelatives)] %>%
c(.[, ID1], .[,ID2]) %>%
.[!(. %in% multiRelatives)]

recovered %>% length

relatedFilt1 <- related[!(ID1 %in% multiRelatives | ID2 %in% multiRelatives)]

relatedSubset <- c(relatedFilt1[,ID1], relatedFilt1[,ID2])

set.seed(4860)

keep <- rbinom(nrow(relatedFilt1), 1,.5)

removeRandom <- relatedFilt1[, remove := ifelse(keep==1, ID1, ID2)] %>%
  .[, remove]


finalSet <- sampQCPheno[!(patIDhda %in% multiRelatives | patIDhda %in% removeRandom)]

finalSet %>% .[,  .(whiteBrit, selfRepEthnicGroup)] %>% table(., useNA = "always")

## Unrelated individuals in other ancestry groups to be defined when ancestry groupings attained.

# Define final discover set and related "replication" set
sampQCPhenoWhiteBrit <- sampQCPheno[patIDhda %in% sampleFile[,ID_1]]

relatedWhiteBrit <- relatedness[ID1 %in% sampQCPhenoWhiteBrit[,patIDhda] & ID2 %in% sampQCPhenoWhiteBrit[,patIDhda]]

relatedWhiteBritIDs <- c(relatedWhiteBrit[,ID1], relatedWhiteBrit[,ID2]) %>%
 unique %>%
 data.table(ID = .)

unrelatedWhiteBritIDs <- sampQCPhenoWhiteBrit[!patIDhda %in% relatedWhiteBritIDs[,ID], patIDhda] %>%
 data.table(ID = .)


fwrite(relatedWhiteBritIDs, file = paste0(outDir,"relatedWhiteBritIDs.txt"))
fwrite(unrelatedWhiteBritIDs, file = paste0(outDir,"unrelatedWhiteBritIDs.txt"))

# non-white british, unrelated set
sampQCPhenoNonWB <- sampQCPheno[is.na(whiteBrit)]

relatedNonWB <- relatedness[ID1 %in% sampQCPhenoNonWB[,patIDhda] & ID2 %in% sampQCPhenoNonWB[,patIDhda]]

relatedWhiteBritIDs <- c(relatedNonWB[,ID1], relatedNonWB[,ID2]) %>%
 unique %>%
 data.table(ID = .)

unrelatedNonWBIDs <- sampQCPhenoNonWB[!patIDhda %in% relatedWhiteBritIDs[,ID], patIDhda] %>%
 data.table(ID = .)


fwrite(unrelatedNonWBIDs, file = paste0(outDir,"unrelatedNonWhiteBritIDs.txt"))
