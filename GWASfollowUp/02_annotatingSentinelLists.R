library(data.table)
library(magrittr)
library(tidyverse)
library(corrplot)
library(patchwork)
library(UpSetR)
library(stringr)



## annotate pixel wise results with FPC hits, and previous hits.

pixelSecondary <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/cojoAnalysis/output/bonfSigLociWithSecondary.csv") %>%
 select(., c("CHR",  "POS",  "SNP_secondary", "POS_secondary", "A1_seconday", 
             "r2_secondary", "BETA_secondary", "SE_secondary", "P_secondary",
             "BETA_secondary_conditional", "SE_secondary_conditional",  
             "P_secondary_conditional")) %>%
  .[, CHR := as.character(CHR)] %>%
  .[CHR=="23", CHR := "X"]

fpcSecondary <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/cojoAnalysis/output/bonfSigFPCLociWithSecondary.csv")%>%
  select(., c("CHR",  "POS",  "SNP_secondary", "POS_secondary", "A1_seconday", 
              "r2_secondary", "BETA_secondary", "SE_secondary", "P_secondary",
              "BETA_secondary_conditional", "SE_secondary_conditional",  
              "P_secondary_conditional")) %>%
  .[, CHR := as.character(CHR)] %>%
  .[CHR=="23", CHR := "X"]


fpcSentinels <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/sentinels/allChr_sentinel_clumpThresh0.001_withOverlap.txt") %>%
  fpcSecondary[., on = c("CHR",  "POS")] %>%
  .[nSNPsLocus >= 5]

sentinels <- lapply(c(1:22, "X"), function(chr) {
  
  chrSent <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/chr",chr,"sentinels_clumpThresh0.001_withOverlap.txt") %>%
    fread() %>%
    .[, CHR := as.character(chr)] %>%
    pixelSecondary[., on =  c("CHR","POS")] %>%
    .[nSNPsLocus >= 5]
  return(chrSent)
}) %>%
  rbindlist()


allPixelsnps <-  c(sentinels[,ID], sentinels[,SNPsInLocus]) %>% 
  strsplit(., ",") %>%
  unlist

allFPCsnps <-  c(fpcSentinels[,ID], fpcSentinels[,SNPsInLocus]) %>% 
  strsplit(., ",") %>%
  unlist
  
rsids <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/misc/rsids_Gao.txt") %>%
  .[, start := NULL]
gao <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/misc/macularThicknessLoci_Gao.txt") %>%
rsids[., on = c("chr", "pos")]

currant1 <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/misc/retinalThicknessLoci_Currant.csv")
currant2 <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/misc/retinalThicknessLoci_Currant_2023.txt")

zekavat <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/misc/Zekavat/genomeWideSigZekavat2024.csv") 

## pixelwise
vep <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/misc/allPixelsSentelsAllVEP_20230501.txt") %>%
  .[, c("#Uploaded_variation", "Consequence", "SYMBOL", "Existing_variation", "DISTANCE", "AF", "PUBMED", "PHENOTYPES")] %>%
  setnames(., "#Uploaded_variation", "ID") %>%
  .[, rsID := sapply(str_split(Existing_variation, pattern = ","), `[`, 1)] %>%
  .[, PUBMED := str_replace_all(PUBMED, ",", ";")] %>%
  .[, PHENOTYPES := str_replace_all(PHENOTYPES, ",", ";")]


freqs <- lapply(c(1:22, "X"), function(chr) {
  fread(paste0("/vast/scratch/users/jackson.v/retThickness/GWAS/annot/chr",chr,"sentinelsFreq.afreq"))
}) %>%
rbindlist %>%
.[, .(ID, REF, ALT, ALT_FREQS)]

snps <- sentinels[, ID]

pixelsAnnotated <- lapply(snps, function(id) {

result <- sentinels[ID %in% id]

allLocusSNPs <-  c(result[,ID], result[,SNPsInLocus]) %>%
  strsplit(., ",") %>%
  unlist

fpcMatch <- fpcSentinels[ID %in% allLocusSNPs]

if(nrow(fpcMatch) > 0) {

  fpcOut <- data.table(fpcSig = "Y",
    topFPC = paste(fpcMatch[,FPC], collapse = ";"),
    allFPCs = paste(fpcMatch[, clumpFPCs], collapse = ";"),
    sameSentinel = ifelse(id %in% fpcMatch[,ID], "Y", "N"),
    FPCsentinel = paste(fpcMatch[,ID], collapse = ";"),
    FPCtopP = paste(fpcMatch[,P], collapse = ";"))

} else {
  fpcOut <- data.table(fpcSig = "N",
  topFPC = NA,
  allFPCs = NA,
  sameSentinel = NA,
  FPCsentinel = NA,
  FPCtopP = NA)

}

gaoMatch <- gao[rsid %in% allLocusSNPs]

if(nrow(gaoMatch) > 0) {

inGao <- data.table(gaoSig = "Y",
    sameSentinelGao = ifelse(id %in% gaoMatch[,rsid], "Y", "N"))

} else {
  inGao <- data.table(gaoSig = "N",
    sameSentinelGao = NA)
}

currantInnerMatch <- currant1[SNP %in% allLocusSNPs]

if(nrow(currantInnerMatch) > 0) {

inCurrantInner <- data.table(currantInnerSig = "Y",
    sameSentinelCurrantInner = ifelse(id %in% currantInnerMatch[,SNP], "Y", "N"))

} else {
inCurrantInner <- data.table(currantInnerSig = "N",
    sameSentinelCurrantInner = NA)
}

currantOuterMatch <- currant2[SNP %in% allLocusSNPs]

if(nrow(currantOuterMatch) > 0) {

inCurrantOuter <- data.table(currantOuterSig = "Y",
    sameSentinelCurrantOuter = ifelse(any(currantOuterMatch[,SNP] %in% id), "Y", "N"))

} else {
inCurrantOuter <- data.table(currantOuterSig = "N",
    sameSentinelCurrantOuter = NA)
}

zekavatMatch <- zekavat[rsid %in% allLocusSNPs]

if(nrow(zekavatMatch) > 0) {
  
  inZekavat <- data.table(zekavatSig = "Y",
                               sameSentinelZekavat = ifelse(any(zekavatMatch[,rsid] %in% id), "Y", "N"))
  
} else {
  inZekavat <- data.table(zekavatSig = "N",
                               sameSentinelZekavat = NA)
}

snpResult <- sentinels[ID==id] %>%
setnames(., "pixel", "topPixel") %>%
  .[, BonferroniSig := ifelse(P < 5E-8/29041, "Y", "N")]

out <- cbind(snpResult, fpcOut, inGao, inCurrantOuter, inCurrantInner, inZekavat) %>%
  as.data.table

if(nrow(fpcOut) > 1){ print(paste(id,"fpc"))}
if(nrow(inGao) > 1){ print(paste(id,"gao"))}
if(nrow(inCurrantInner) > 1){ print(paste(id,"currantInner"))}
if(nrow(inCurrantOuter) > 1){ print(paste(id,"currantOuter"))}
if(nrow(inZekavat) > 1){ print(paste(id,"zekavatOuter"))}

return(out)

}) %>%
rbindlist %>%
.[, novel := ifelse(gaoSig=="N" & currantInnerSig=="N" & currantOuterSig=="N" & zekavatSig=="N", "Y", "N")] %>% 
.[vep, on = "ID"] %>%
.[freqs, on = "ID"] %>%
.[!is.na(POS)] %>%  
.[, EffAlleleFreq := ifelse(A1==ALT, ALT_FREQS, (1-ALT_FREQS))] %>%
.[, EffAllele := A1] %>%
.[, nonEffAllele := ifelse(A1==ALT, REF, ALT)] %>%
.[, .(CHR, POS, ID, rsID, EffAllele, nonEffAllele, EffAlleleFreq, BETA, SE, P, BonferroniSig, topPixel,
nPixelsLocus, nSNPsLocus, fpcSig, topFPC, allFPCs, sameSentinel, FPCsentinel,
FPCtopP, gaoSig, sameSentinelGao, currantOuterSig, sameSentinelCurrantOuter,
currantInnerSig, sameSentinelCurrantInner, zekavatSig, 
sameSentinelZekavat, novel,  Consequence, SYMBOL,
novel,Consequence, SYMBOL, DISTANCE, AF, PUBMED, PHENOTYPES, SNP_secondary,
POS_secondary, A1_seconday, r2_secondary, BETA_secondary, SE_secondary,
P_secondary, BETA_secondary_conditional, SE_secondary_conditional,  
P_secondary_conditional)] %>%
  setnames(., "AF", "AF1000G") %>%
  unique





## FPC results

vepFPC <-  fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/misc/FPCsSentelsAllVEP_20230501.txt") %>%
  .[, c("#Uploaded_variation", "Consequence", "SYMBOL", "Existing_variation", "DISTANCE", "AF", "PUBMED", "PHENOTYPES")] %>%
  setnames(., "#Uploaded_variation", "ID") %>%
  .[, rsID := sapply(str_split(Existing_variation, pattern = ","), `[`, 1)] %>%
  .[, PUBMED := str_replace_all(PUBMED, ",", ";")] %>%
  .[, PHENOTYPES := str_replace_all(PHENOTYPES, ",", ";")]


snps <- fpcSentinels[, ID]

fpcsAnnotated <- lapply(snps, function(id) {

result <- fpcSentinels[ID==id]

allLocusSNPs <-  c(result[,ID], result[,SNPsInLocus]) %>% 
  strsplit(., ",") %>%
  unlist

pixelMatch <- sentinels[ID %in% allLocusSNPs]

if(nrow(pixelMatch) > 0) {

  pixelOut <- data.table(pixelwiseSig = "Y",
    toppixel = paste(pixelMatch[,pixel], collapse = ";"),
    nPixels = paste(pixelMatch[, nPixelsLocus], collapse = ";"),
    sameSentinel = ifelse(id %in% pixelMatch[,ID], "Y", "N"),
    pixelSentinel = paste(pixelMatch[,ID], collapse = ";"),
    pixelTopP = paste(pixelMatch[,P], collapse = ";"))

} else {
  pixelOut <- data.table(pixelwiseSig = "N",
  toppixel = NA,
  nPixels = NA,
  sameSentinel = NA,
  pixelSentinel = NA,
  pixelTopP = NA)

}

gaoMatch <- gao[rsid %in% allLocusSNPs]

if(nrow(gaoMatch) > 0) {

inGao <- data.table(gaoSig = "Y",
    sameSentinelGao = ifelse(id %in% gaoMatch[,rsid], "Y", "N"))

} else {
  inGao <- data.table(gaoSig = "N",
    sameSentinelGao = NA)
}

currantInnerMatch <- currant1[SNP %in% allLocusSNPs]

if(nrow(currantInnerMatch) > 0) {

inCurrantInner <- data.table(currantInnerSig = "Y",
    sameSentinelCurrantInner = ifelse(id %in% currantInnerMatch[,SNP], "Y", "N"))

} else {
inCurrantInner <- data.table(currantInnerSig = "N",
    sameSentinelCurrantInner = NA)
}

currantOuterMatch <- currant2[SNP %in% allLocusSNPs]

if(nrow(currantOuterMatch) > 0) {

inCurrantOuter <- data.table(currantOuterSig = "Y",
    sameSentinelCurrantOuter = ifelse(any(currantOuterMatch[,SNP] %in% id), "Y", "N"))

} else {
inCurrantOuter <- data.table(currantOuterSig = "N",
    sameSentinelCurrantOuter = NA)
}

zekavatMatch <- zekavat[rsid %in% allLocusSNPs]

if(nrow(zekavatMatch) > 0) {
  
  inZekavat <- data.table(zekavatSig = "Y",
                          sameSentinelZekavat = ifelse(any(zekavatMatch[,rsid] %in% id), "Y", "N"))
  
} else {
  inZekavat <- data.table(zekavatSig = "N",
                          sameSentinelZekavat = NA)
}

snpResult <- fpcSentinels[ID==id] %>%
  setnames(., c("FPC", "clumpFPCs"), c("topFPC", "sigFPCs")) %>%
  .[, BonferroniSig := ifelse(P < 5E-8/6, "Y", "N")]

out <- cbind(snpResult, pixelOut, inGao, inCurrantOuter, inCurrantInner, inZekavat) %>%
  as.data.table

if(nrow(pixelOut) > 1){ print(paste(id,"pixels"))}
if(nrow(inGao) > 1){ print(paste(id,"gao"))}
if(nrow(inCurrantInner) > 1){ print(paste(id,"currantInner"))}
if(nrow(inCurrantOuter) > 1){ print(paste(id,"currantOuter"))}
if(nrow(inZekavat) > 1){ print(paste(id,"currantOuter"))}

return(out)

}) %>%
rbindlist %>%
.[, novel := ifelse(gaoSig=="N" & currantInnerSig=="N" & currantOuterSig=="N" & zekavatSig=="N", "Y", "N")] %>% 
.[vepFPC, on = "ID"] %>%
  .[!is.na(POS)] %>%  
  .[, EffAlleleFreq := A1_FREQ] %>%
.[, EffAllele := A1] %>%
.[, nonEffAllele := ifelse(A1==ALT, REF, ALT)] %>%
.[, .(CHR, POS, ID, rsID, EffAllele, nonEffAllele, EffAlleleFreq, BETA, SE, P, BonferroniSig, topFPC,
sigFPCs, nSNPsLocus, pixelwiseSig, toppixel, nPixels, sameSentinel, pixelSentinel,
pixelTopP, gaoSig, sameSentinelGao, currantOuterSig, sameSentinelCurrantOuter,
currantInnerSig, sameSentinelCurrantInner, zekavatSig, sameSentinelZekavat, 
novel,  Consequence, SYMBOL, DISTANCE, AF, PUBMED, PHENOTYPES, SNP_secondary,
POS_secondary, A1_seconday, r2_secondary, BETA_secondary, SE_secondary,
P_secondary, BETA_secondary_conditional, SE_secondary_conditional,  
P_secondary_conditional)] %>%
setnames(., "AF", "AF1000G") %>%
  unique



fwrite(pixelsAnnotated, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/annotatedSentinelsPixelwise.txt", sep = "\t")
fwrite(pixelsAnnotated[novel=="Y"], file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/annotatedSentinelsPixelwiseNovelOnly.txt", sep = "\t")
fwrite(pixelsAnnotated[BonferroniSig=="Y", .(rsID)], file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/rsIDsSentinelsPixelwiseBonferroniSig.txt", sep = "\t")

fwrite(fpcsAnnotated, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/annotatedSentinelsFPCs.txt", sep = "\t")
fwrite(fpcsAnnotated[novel=="Y"], file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/annotatedSentinelsFPCsNovelOnly.txt", sep = "\t")
fwrite(fpcsAnnotated[BonferroniSig=="Y", .(rsID)], file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/rsIDsSentinelsFPCsBonferroniSig.txt", sep = "\t")


clusters <- fread("/wehisan/bioinf/lab_bahlo/projects/mactel/PRJ2022002/analysis_versions/version001/Post-GWAS/results/cluster_intepretation.csv", fill =T)

# file fucked
# Function to extract correct info
extractCols <- function(row) {
  final_index <- max(which(!is.na(row) & row != ""))
  return(c(row[3], row[final_index]))
}

# Apply the function to each row of the data frame
clusterAssignments <- t(apply(clusters, 1, extractCols)) %>%
as.data.table %>%
setnames(., c("ID", "clusterAssignment")) 

## output loci summary for supp table
pixelsSuppTab <- merge.data.table(pixelsAnnotated, clusterAssignments, on = c("ID"), all.x = T) %>%
                  .[, .(CHR, POS, ID, rsID, EffAllele, nonEffAllele, EffAlleleFreq, BETA, SE, P, BonferroniSig, topPixel,
                                   nPixelsLocus, nSNPsLocus, fpcSig, topFPC, allFPCs, sameSentinel, FPCsentinel,
                                   FPCtopP, novel, SNP_secondary,  POS_secondary, A1_seconday, r2_secondary, BETA_secondary, 
                                   SE_secondary, P_secondary, BETA_secondary_conditional, SE_secondary_conditional,  
                                   P_secondary_conditional, clusterAssignment)]

fpcsSuppTab <- fpcsAnnotated[, .(CHR, POS, ID, rsID, EffAllele, nonEffAllele, EffAlleleFreq, BETA, SE, P, BonferroniSig, topFPC,
                                      sigFPCs, nSNPsLocus, pixelwiseSig, toppixel, nPixels, sameSentinel, pixelSentinel,
                                      pixelTopP, novel, SNP_secondary,  POS_secondary, A1_seconday, r2_secondary, BETA_secondary, 
                                      SE_secondary, P_secondary, BETA_secondary_conditional, SE_secondary_conditional,  
                                      P_secondary_conditional)]


fwrite(pixelsSuppTab , file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/annotatedLociPixelwise.txt", sep = "\t")
fwrite(fpcsSuppTab, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/annotatedLociFPCs.txt", sep = "\t")


## output overlap tables for supplement
pixOverlap <- pixelsAnnotated %>%
  .[BonferroniSig == "Y"] %>%
  .[, .(CHR, POS, ID, rsID, EffAllele, nonEffAllele, EffAlleleFreq, BETA, SE, P, BonferroniSig, topPixel,
                    nSNPsLocus, fpcSig, topFPC, sameSentinel, gaoSig, sameSentinelGao, currantOuterSig, sameSentinelCurrantOuter,
                    currantInnerSig, sameSentinelCurrantInner, zekavatSig, sameSentinelZekavat, 
                    novel)] %>%
  unique 


fpcOverlap <- fpcsAnnotated %>%
  .[BonferroniSig == "Y"] %>%
  .[, .(CHR, POS, ID, rsID, EffAllele, nonEffAllele, EffAlleleFreq, BETA, SE, P, BonferroniSig, topFPC,
                                  nSNPsLocus, pixelwiseSig, toppixel, sameSentinel, 
                                  gaoSig, sameSentinelGao, currantOuterSig, sameSentinelCurrantOuter,
                                  currantInnerSig, sameSentinelCurrantInner, zekavatSig, sameSentinelZekavat, 
                                  novel)] %>%
  unique 



fwrite(pixOverlap , file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/overlapPreviousLociPixelwise.txt", sep = "\t")

fwrite(fpcOverlap , file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations2024-05/overlapPreviousLociFPCs.txt", sep = "\t")

