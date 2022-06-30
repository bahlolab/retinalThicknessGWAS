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
scansDTlinked <- linkage[scansDT, on =  "patID"]
# scansDTlinked <- linkage[scansDT, on = c("patIDhda" = "patID")]

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

#write covariate file- both with single, and double IDs
fwrite(covsOut, file = paste0(outDir,"covariates_doubleIDs.txt"), sep = "\t", na = "NA", quote = F)
fwrite(covsOut[,!"FID"], file = paste0(outDir,"covariates_singleIDs.txt"), sep = "\t", na = "NA", quote = F)

phenoOut <- scansDTlinked[, !c("patID", "visit", "eye", "sex", "age", "device", "meanRefErr")] %>%
  cbind(data.table(FID = .[,patIDhda], IID = .[,patIDhda]), .) %>%
  .[, patIDhda := NULL]


pixels <- data.table(pixel = names(phenoOut)[!names(phenoOut) %in% c("FID", "IID")]) %>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert=TRUE)]

fwrite(pixels, file = paste0(outDir,"pixels.txt"), sep = "\t", na = "NA", quote = F, col.names = F)

slices <- pixels[,y] %>% unique

sapply(slices, function(idx) {

  cols <- c("FID", "IID", pixels[y==idx, pixel])
  out <- phenoOut[, ..cols]

  fwrite(out, file = paste0(outDir,"phenotypesSlice",idx,"_doubleIDs.txt"), sep = "\t", na = "NA", quote = F)
  fwrite(out[,!"FID"], file = paste0(outDir,"phenotypesSlice",idx,"_singleIDs.txt"), sep = "\t", na = "NA", quote = F)

})
