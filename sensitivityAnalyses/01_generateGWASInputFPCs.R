## Run using R/4.1.2


library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)
library(DescTools)
library(fastDummies)

dataDir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/sensitivityAnalysesGWAS/rawData/"
outDir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/sensitivityAnalysesGWAS/output/"

scansDT  <- fread(paste0(dataDir,"scansUnadjustedFinal.csv"), select = c("patID", "visit", "eye", "sex", "age", "device", "meanRefErr"))
fpcDT <- fread(paste0(dataDir,"fPCscores_noExclusions.csv"))


linkage <- fread(paste0(dataDir,"idLinkage.txt"))
scansDTlinked <- linkage[fpcDT, on = "patID"] %>%
  .[scansDT, on =  "patID"] %>%
  .[!is.na(patIDhda)]




#####################################################
## identify individuals who underwent surgery.
surgery <- fread(paste0(dataDir,"ukb11226.tab"),
                  select=c("f.eid", "f.5181.0.0", "f.5181.1.0", "f.5324.0.0", "f.5324.1.0", "f.5325.0.0", "f.5325.1.0", "f.5326.0.0", "f.5326.1.0", "f.5327.0.0", "f.5327.1.0", "f.5328.0.0", "f.5328.1.0"), quote="") %>%
  .[, glaucomaSurgery := as.integer(rowSums(.SD, na.rm = TRUE) > 0), .SDcols = c("f.5326.0.0", "f.5326.1.0", "f.5327.0.0", "f.5327.1.0")] %>%
  .[, glaucomaSurgeryDefinite := as.integer(rowSums(.SD, na.rm = TRUE) > 0 & rowSums(.SD, na.rm = TRUE) <5), .SDcols = c("f.5326.0.0", "f.5326.1.0", "f.5327.0.0", "f.5327.1.0")] %>%
  .[, .( f.eid, glaucomaSurgery, glaucomaSurgeryDefinite)] %>%
  .[linkage, on = c("f.eid" = "patID")]


## make smoking pheno file
smoke<- data.table(f.eid = scansDTlinked[,patID])

smokeInstances <- lapply(c(0,1), function(instance) {

instanceDT <- lapply(c("current","previous", "never"), function(status) {
  
  dt<- paste0(dataDir,"/",status,"Smoke_",instance,".csv") %>%
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
.[smoke[, .(f.eid, smokingFactor_0, smokingFactor_1)], on = "f.eid"]

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
  .[, c("patIDhda", "f.50.0.0", "f.50.1.0", "f.eid", "visit", "device", "smokingFactor_0", "smokingFactor_1") := NULL] %>%
  .[!is.na(IID)]


phenoFull <- scansDTlinked[, !c("patID", "visit", "eye", "sex", "age", "device", "meanRefErr")] %>%
  cbind(data.table(FID = .[,patIDhda], IID = .[,patIDhda]), .) %>%
  .[, patIDhda := NULL]


# fwrite(covsFull, file = paste0(outDir,"covariatesFullCohort.txt"), sep = "\t", na = "NA", quote = F)
fwrite(phenoFull, file = paste0(outDir,"FPCphenotypesFullCohort.txt"), sep = "\t", na = "NA", quote = F)


anc <- "EUR"
  
  ids <- paste0(dataDir,"sampleList_singleIDs_",anc,".txt") %>%
    fread
  
  print(paste(nrow(ids), "individuals with", anc, "ancestry."))
  print(paste(nrow(covsFull[IID %in% ids[,IID]]), "individuals with covariate data."))
  print(paste(nrow(phenoFull[IID %in% ids[,IID]]), "individuals with phenotype data."))
  

## Sensitivity analyses:
## 1. add smoking as a covariate
## 2. exclude individuals who underwent glaucoma surgery

## 1. smoking as a covariate
smokingCovs <- covsFull[IID %in% ids[,IID]] %>%
  .[!is.na(smokingStatus)]

print(paste(nrow(smokingCovs), "individuals with smoking data."))

smokingPheno <- phenoFull[IID %in% smokingCovs[,IID]]


# fwrite(smokingCovs, file = paste0(outDir,"covariates_doubleIDs_smoking_",anc,".txt"), sep = "\t", na = "NA", quote = F)
fwrite(smokingPheno, file = paste0(outDir,"FPCphenotypes_doubleIDs_smoking_",anc,".txt"), sep = "\t", na = "NA", quote = F)


rm(smokingCovs, smokingPheno)

## 2. exclude individuals who underwent glaucoma surgery
hadSurgeryDefinite <- surgery[glaucomaSurgeryDefinite == 1, patIDhda]

noSurgeryCovs <- covsFull[IID %in% ids[,IID]] %>%
  .[!IID %in% hadSurgeryDefinite] %>%
  .[, smokingStatus := NULL]

print(paste(nrow(noSurgeryCovs), "individuals without glaucoma surgery."))

noSurgeryPheno <- phenoFull[IID %in% noSurgeryCovs[,IID]]

# fwrite(noSurgeryCovs, file = paste0(outDir,"covariates_doubleIDs_noSurgery_",anc,".txt"), sep = "\t", na = "NA", quote = F)
fwrite(noSurgeryPheno, file = paste0(outDir,"FPCphenotypes_doubleIDs_noSurgery_",anc,".txt"), sep = "\t", na = "NA", quote = F)
