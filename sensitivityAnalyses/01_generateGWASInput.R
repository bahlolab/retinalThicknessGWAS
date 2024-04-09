## Run using R/4.1.2


library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)
library(DescTools)
library(fastDummies)

dataDir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/rawData/"
outDir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/"

scansDT  <-fread(paste0(dataDir,"scansUnadjustedFinal.csv"))


linkage <- fread(paste0(dataDir,"idLinkage.txt"))
scansDTlinked <- linkage[scansDT, on =  "patID"] %>%
  .[!is.na(patIDhda)]
# scansDTlinked <- linkage[scansDT, on = c("patIDhda" = "patID")]




#####################################################
## identify individuals who underwent surgery.
surgery <- fread("/stornext/Bioinf/data/lab_bahlo/projects/misc/UKBiobank/data/app28541/phenoData/ukb11226.tab",
                  select=c("f.eid", "f.5181.0.0", "f.5181.1.0", "f.5324.0.0", "f.5324.1.0", "f.5325.0.0", "f.5325.1.0", "f.5326.0.0", "f.5326.1.0", "f.5327.0.0", "f.5327.1.0", "f.5328.0.0", "f.5328.1.0"), quote="") %>%
  .[, glaucomaSurgery := as.integer(rowSums(.SD, na.rm = TRUE) > 0), .SDcols = c("f.5326.0.0", "f.5326.1.0", "f.5327.0.0", "f.5327.1.0")] %>%
  .[, glaucomaSurgeryDefinite := as.integer(rowSums(.SD, na.rm = TRUE) > 0 & rowSums(.SD, na.rm = TRUE) <5), .SDcols = c("f.5326.0.0", "f.5326.1.0", "f.5327.0.0", "f.5327.1.0")] %>%
  .[, .( f.eid, glaucomaSurgery, glaucomaSurgeryDefinite)]

hadSurgeryDefinite <- surgery[glaucomaSurgeryDefinite == 1, f.eid]

## make smoking pheno file
smokeDir <- "/stornext/Bioinf/data/lab_bahlo/projects/misc/UKBiobank/data/app28541/phenoData/smokingStatusIDs"
smoke<- data.table(f.eid = scansDTlinked[,patID])

smokeInstances <- lapply(c(0,1), function(instance) {

instanceDT <- lapply(c("current","previous", "never"), function(status) {
  
  dt<- paste0(smokeDir,"/",status,"Smoke_",instance,".csv") %>%
  fread %>%
  setnames(., c("f.eid")) %>%
  .[, status := status]

return(dt)  
}) %>%
rbindlist %>%
setnames(., "status", paste0("smokingStatus_",instance))

})

## merge smokeInstances with smoke, on f.eid
## recode to factor for each instance to never=0, previous=1, current=2
smoke <-    merge(smoke, smokeInstances[[1]], by = "f.eid", all.x = TRUE) %>%
  merge(., smokeInstances[[2]], by = "f.eid", all.x = TRUE) %>%
  .[, smokingFactor_0 := case_when(smokingStatus_0 == "never" ~ 0,
                                 smokingStatus_0 == "previous" ~ 1,
                                 smokingStatus_0 == "current" ~ 2 )] %>%
  .[, smokingFactor_1 := case_when(smokingStatus_1 == "never" ~ 0,
                                  smokingStatus_1 == "previous" ~ 1,
                                  smokingStatus_1 == "current" ~ 2 )] 








file <- paste0(dataDir,"ukb41258.tab")
names <- fread(file, nrows = 0) %>%
  names

getIdx <- function(code) {
    names %like any% code %>% which
  }

cols <- getIdx(c("f.eid", "f.50.0.0", "f.50.1.0"))

covs <- fread(file, select = cols) %>%
.[smoke[, .(f.eid, smokingFactor_0, smokingFactor_1)], on = f.eid]

pcs <- paste0(dataDir,"ukb32825.tab") %>%
  fread(., select = c(1, 1145:1154)) %>%
  setnames(. ,c("patIDhda", paste0("PC",c(1:10))))

covsFull <- scansDTlinked[, c("patIDhda", "patID", "visit", "eye", "sex", "age", "device", "meanRefErr")] %>%
  covs[., on = c(f.eid = "patID")] %>%
  pcs[., on = "patIDhda"] %>%
  .[, standHeight := case_when(visit == 0 ~ f.50.0.0,
                      visit == 1 ~ f.50.1.0)] %>%
  .[, ageSq := age^2] %>%
  .[, smokingStatus := case_when(visit == 0 ~ smokingFactor_0,
                                visit == 1 ~ smokingFactor_1)] %>%
  # fastDummies::dummy_cols(., select_columns = c("eye", "device")) %>%
  .[, deviceCat := paste0("dev",device)] %>%
  cbind(data.table(FID = .[,patIDhda], IID = .[,patIDhda]), .) %>%
  .[, c("patIDhda", "f.50.0.0", "f.50.1.0", "f.eid", "visit", "device") := NULL] %>%
  .[!is.na(IID)]


phenoFull <- scansDTlinked[, !c("patID", "visit", "eye", "sex", "age", "device", "meanRefErr")] %>%
  cbind(data.table(FID = .[,patIDhda], IID = .[,patIDhda]), .) %>%
  .[, patIDhda := NULL]


pixels <- data.table(pixel = names(phenoOut)[!names(phenoOut) %in% c("FID", "IID")]) %>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert=TRUE)]

# fwrite(pixels, file = paste0(outDir,"pixels.txt"), sep = "\t", na = "NA", quote = F, col.names = F)

slices <- pixels[,y] %>% unique

anc <- "EUR"
  
  ids <- paste0(dataDir,"sampleList_singleIDs_",anc,".txt") %>%
    fread
  
  print(paste(nrow(ids), "individuals with", anc, "ancestry."))
  print(paste(nrow(covsOut[IID %in% ids[,IID]]), "individuals with covariate data."))
  print(paste(nrow(phenoOut[IID %in% ids[,IID]]), "individuals with phenotype data."))
  

## Sensitivity analyses:
## 1. add smoking as a covariate
## 2. exclude individuals who underwent glaucoma surgery

#write covariate file- both with single, and double IDs
fwrite(covsOut[IID %in% ids[,IID]], file = paste0(outDir,"covariates_doubleIDs_",anc,".txt"), sep = "\t", na = "NA", quote = F)




sapply(slices, function(idx) {

  cols <- c("FID", "IID", pixels[y==idx, pixel])
  out <- phenoOut[, ..cols]

  fwrite(out[IID %in% ids[,IID]], file = paste0(outDir,"phenotypesSlice",idx,"_doubleIDs_",anc,".txt"), sep = "\t", na = "NA", quote = F)

})

})
