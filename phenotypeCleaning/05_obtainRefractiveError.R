## Run using R/4.1.2


## Get refractive errors
## these will be used for adjusted phnotypes further down the line.

library(data.table)
library(magrittr)
library(dplyr)
library(parallel)


file <- "/wehisan/bioinf/lab_bahlo/projects/misc/UKBiobank/data/app28541/phenoData/ukb41258.tab"
names <- fread(file, nrows = 0) %>%
  names

getIdx <- function(code) {
    names %like% code %>% which
  }


getBest <- function(best, measures) {

  bestIdx <- getIdx(best)
  measuresIdx <- getIdx(measures)

  dt <- fread(file, select = measuresIdx)
  ref <- fread(file, select = bestIdx) %>% unlist

  vec <- c(1:length(ref))

  bestMeasure <- lapply(vec, function(x) {

    r <- ref[x]

    if(!is.na(r)) {

      bestCol <- paste(measures, r, sep = ".")
      meas <-  dt[x, ..bestCol] %>% unlist

    } else {

      meas <- NA

    }
    return(meas)
  })  %>% unlist
}


## obtain best sperical and cyndrilical thing for each eye

cylindL0 <- getBest(best = "f.5276.0", measures = "f.5086.0")
cylindL1 <- getBest(best = "f.5276.1", measures = "f.5086.1")
cylindR0 <- getBest(best = "f.5221.0", measures = "f.5087.0")
cylindR1 <- getBest(best = "f.5221.1", measures = "f.5087.1")

sphereL0 <- getBest(best = "f.5276.0", measures = "f.5085.0")
sphereL1 <- getBest(best = "f.5276.1", measures = "f.5085.1")
sphereR0 <- getBest(best = "f.5221.0", measures = "f.5084.0")
sphereR1 <- getBest(best = "f.5221.1", measures = "f.5084.1")


## table of refractive errors for each eye/instances
ids <- fread(file, select = 1) %>% unlist
refError <- data.table(id = ids,
                      L0 = sphereL0 + 0.5 * cylindL0,
                      L1 = sphereL1 + 0.5 * cylindL1,
                      R0 = sphereR0 + 0.5 * cylindR0,
                      R1 = sphereR1 + 0.5 * cylindR1 )

fwrite(refError, "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/refractiveError.csv", sep = ",")
