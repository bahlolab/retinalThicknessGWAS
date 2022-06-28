
## Run using R/4.1.2

library(data.table)
library(magrittr)
library(dplyr)
library(parallel)
library(rlist)
library(lubridate)
library(ggplot2)
library(plotly)



load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/gridMask.RData")

sampsIDs <- fread("/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/data/scanIDs.txt",
                col.names = c("id", "patientID"))

uniqueSamps <- sampsIDs[,id] %>% unique

summlist <- function(x) {

  which(uniqueSamps==x) %>% paste("Sample no.", .) %>% print

  scans <- sampsIDs[id == x, patientID]

  sampSumm <- lapply(scans, function(y) {

    # print(y)
    file <- paste0("/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/data/rawPerScan/",x,"/",y,"_scan.csv")

    if (file.exists(file))  {
    dt <- fread(file, col.names = c("scan", "eye", "slice_index", "fovea_index", 0:255)) %>%
      .[, c("id", "laterality", "visit", "measure") := tstrsplit(scan, "_")]

    fullScanMat <- dt[, c("slice_index", 0:255)] %>%
      as.matrix(., rownames="slice_index")


    # set -1 to NA
    fullScanMat[fullScanMat == -1] <- NA

    # trim
    fullScanMat <- replace(fullScanMat, missing0.1Final == 1, NA)
    pointsRemoved <- missing0.1Final  %>% sum

    # summarise missingness for each scan
    naCount <- fullScanMat %>% is.na %>% sum %>% subtract(pointsRemoved)

    ## calculate scan-wise summaries
    myMean <- mean(fullScanMat, na.rm = T)
    mySD <- sd(fullScanMat, na.rm = T)

    myMedian <- median(fullScanMat, na.rm = T)


    # data table to retun
    summ <- dt[, c("scan", "eye", "id", "laterality", "visit", "measure")] %>%
        unique %>%
        .[, meanDepth := myMean] %>%
        .[, medianDepth := myMedian] %>%
        .[, minDepth := min(fullScanMat, na.rm = T)] %>%
        .[, maxDepth := max(fullScanMat, na.rm = T)] %>%
        .[, sdDepth := mySD] %>%
        .[, totalMissing := naCount]

    return(summ)
  }}) %>%
  rbindlist
}

scansSumm <- mclapply(uniqueSamps, summlist, mc.cores = 8)

save(scansSumm, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/scansSummariesTrimmed20220503.RData")

fields <- c("meanDepth", "medianDepth", "minDepth",  "maxDepth", "sdDepth", "totalMissing")


lapply(fields, function(f) {

  plot <- scansSumm %>%
    rbindlist %>%
    select(.,  c("scan", all_of(f)))

  ggplot(plot, aes_string(x = f)) +
    geom_histogram() +
    ggtitle(f) +
    scale_y_continuous(trans='log2')

    ggsave(file = paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/plots/trimmedGrid_summaries_",f,"_log.png"))

    ggplot(plot, aes_string(x = f)) +
      geom_histogram() +
      ggtitle(f)

    ggsave(file = paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/plots/trimmedGrid_summaries_",f,".png"))

} )


## Function to select "best" scan id(s).
selectScan <- function(x) {

  samp <- x[,id] %>% unique
  which(uniqueSamps==samp) %>% paste("Sample no.", .) %>% print

  # first remove scans with > ?? datapoints missing
  x <- x[totalMissing < 29041*0.1]

  if(nrow(x) > 0) {

  scans <- vector()

  # counts for each eye / instance
  nL0 <- x[laterality==21011 & visit==0] %>% nrow
  nL1 <- x[laterality==21011 & visit==1] %>% nrow
  nR0 <- x[laterality==21013 & visit==0] %>% nrow
  nR1 <- x[laterality==21013 & visit==1] %>% nrow

  # refractive error available for each eye / instance
  refErrL0 <- refError[id == samp , !is.na(L0)]
  refErrL1 <- refError[id == samp , !is.na(L1)]
  refErrR0 <- refError[id == samp , !is.na(R0)]
  refErrR1 <- refError[id == samp , !is.na(R1)]

  # select best scan per eye, per instance (based on missingness), if more than one
  selectL0 <- selectL1 <- selectR0 <- selectR1 <-vector()

  if(nL0>=1 & refErrL0) {
    selectL0 <- x[laterality==21011 & visit==0] %>%
      .[totalMissing==min(totalMissing)] %>%
      .[, scan] %>%
      sample(., 1) # if ties, randomly select one
  }

  if(nR0>=1 & refErrR0) {
    selectR0 <- x[laterality==21013 & visit==0] %>%
      .[totalMissing==min(totalMissing)] %>%
      .[, scan] %>%
      sample(., 1)
  }

  if(nL1>=1 & refErrL1) {
    selectL1 <- x[laterality==21011 & visit==1] %>%
      .[totalMissing==min(totalMissing)] %>%
      .[, scan] %>%
      sample(., 1)
  }

  if(nR1>=1 & refErrR1) {
    selectR1 <- x[laterality==21013 & visit==1] %>%
      .[totalMissing==min(totalMissing)] %>%
      .[, scan] %>%
      sample(., 1)
  }

  # both eyes measured at each instance?

  both0 <- ifelse(nL0>=1 & nR0>=1 & refErrL0 & refErrR0, 1, 0)
  both1 <- ifelse(nL1>=1 & nR1>=1 & refErrL1 & refErrR1, 1, 0)

  ## both eyes, at one instance
  if(both0==1 & both1==0) { scans <- c(selectL0, selectR0)}
  if(both0==0 & both1==1) { scans <- c(selectL1, selectR1)}

  ## both eyes,both instances
  if(both0==1 & both1==1) {

    missing0 <- x[scan %in% c(selectL0, selectR0)] %>%
      .[, totalMissing] %>% sum

    missing1 <- x[scan %in% c(selectL1, selectR1)] %>%
      .[, totalMissing] %>% sum

    scans <- ifelse(missing0 <= missing1, c(selectL0, selectR0),  c(selectL1, selectR1))

  }

  ## both eyes, neither instance
  if(both0==0 & both1==0 &
          length(c(selectL0, selectR0, selectL1, selectR1)) > 0) {
    scans <- x %>%
      .[scan %in% c(selectL0, selectR0, selectL1, selectR1)] %>%
      .[totalMissing==min(totalMissing)] %>%
      .[, scan] %>%
      sample(., 1) # if ties, randomly select one
  }


  return(scans)
}
}


## Get available scans, by eye and instance
scansAvail <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/scansAvailable20220427.txt")

## look at combinations to check all covered by function...
allCombinations <- scansAvail[, combo := paste0(L_0,R_0,L_1,R_1)] %>% .[,combo] %>% table

## read in refractive error
refError <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/refractiveError.csv")

## select best scan
set.seed(63259)
bestScanIdx <- mclapply(scansSumm, selectScan, mc.cores = 8)

save(bestScanIdx, file="/stornext/Bioinf/data/lab_bahlo/projects/misc/retinalThickness/phenoExploratory/bestScansIndex20220504.RData")
# x <- scansSumm[[1]]

lapply(bestScanIdx, length) %>% unlist %>% table
