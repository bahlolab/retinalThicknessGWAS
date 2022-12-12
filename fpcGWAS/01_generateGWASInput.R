## Run using R/4.1.2


library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)
library(DescTools)
library(fastDummies)

dataDir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWAS/rawData/"
outDir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWAS/output/"

scansDT  <- fread(paste0(dataDir,"scansUnadjustedFinalfPCexclusions.csv"), select = c("patID", "visit", "eye", "sex", "age", "device", "meanRefErr"))
fpcDT <- fread(paste0(dataDir,"fpcPhenotypes.txt"))


linkage <- fread(paste0(dataDir,"idLinkage.txt"))
scansDTlinked <- linkage[fpcDT, on = "patID"] %>%
  .[scansDT, on =  "patID"] %>%
  .[!is.na(patIDhda)]

file <- paste0(dataDir,"ukb41258.tab")
names <- fread(file, nrows = 0) %>%
  names

getIdx <- function(code) {
    names %like any% code %>% which
  }

cols <- getIdx(c("f.eid", "f.50.0.0", "f.50.1.0"))

covs <- fread(file, select = cols)

pcs <- paste0(dataDir,"ukb32825.tab") %>%
  fread(., select = c(1, 1145:1154)) %>%
  setnames(. ,c("patIDhda", paste0("PC",c(1:10))))

covsOut <- scansDTlinked[, c("patIDhda", "patID", "visit", "eye", "sex", "age", "device", "meanRefErr")] %>%
  covs[., on = c(f.eid = "patID")] %>%
  pcs[., on = "patIDhda"] %>%
  .[, standHeight := case_when(visit == 0 ~ f.50.0.0,
                      visit == 1 ~ f.50.1.0)] %>%
  .[, ageSq := age^2] %>%
  # fastDummies::dummy_cols(., select_columns = c("eye", "device")) %>%
  .[, deviceCat := paste0("dev",device)] %>%
  cbind(data.table(FID = .[,patIDhda], IID = .[,patIDhda]), .) %>%
  .[, c("patIDhda", "f.50.0.0", "f.50.1.0", "f.eid", "visit", "device") := NULL] %>%
  .[!is.na(IID)]

phenoOut <- scansDTlinked[, !c("patID", "visit", "eye", "sex", "age", "device", "meanRefErr")] %>%
  cbind(data.table(FID = .[,patIDhda], IID = .[,patIDhda]), .) %>%
  .[, patIDhda := NULL]


lapply(c("EUR", "CSA", "AFR"), function(anc) {
 
  ids <- paste0(dataDir,"sampleList_singleIDs_",anc,".txt") %>%
    fread
  
  print(paste(nrow(ids), "individuals with", anc, "ancestry."))
  print(paste(nrow(covsOut[IID %in% ids[,IID]]), "individuals with covariate data."))
  print(paste(nrow(phenoOut[IID %in% ids[,IID]]), "individuals with phentype data."))
  
   #write covariate file- both with single, and double IDs
  fwrite(covsOut[IID %in% ids[,IID]], file = paste0(outDir,"covariates_doubleIDs_",anc,".txt"), sep = "\t", na = "NA", quote = F)

  
  fwrite(phenoOut[IID %in% ids[,IID]], file = paste0(outDir,"FPCphenotypes_doubleIDs_",anc,".txt"), sep = "\t", na = "NA", quote = F)


})



