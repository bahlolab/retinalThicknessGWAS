
## Run using R/4.1.2

## Generate summaries on raw phenotype data.
## Basically a sense-check step, to confirm all looks as expected...

library(data.table)
library(magrittr)
library(dplyr)
library(parallel)
library(rlist)
library(lubridate)
library(ggplot2)
library(plotly)


## Read in list of healthy individuals -  use list from Yue/Yuka
phenoDir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/data/phenotypes/rawData/"

healthyFinal <- fread(paste0(phenoDir,"normal_patients.txt")) %>%
  .[,V1]


## list of scans to remove as outliers
moorFieldsExclude <- fread(paste0(phenoDir,"exclusion_manual_Moorfiled.csv"))
outliers <- fread(paste0(phenoDir,"outliers_all.csv"), header=F)

sliceDir <- "/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/data/rawData/"

## read in middle foveal scan to get list of ids
scanIdx <- 63

dt <- paste0(sliceDir,"height_",scanIdx,".csv") %>%
    fread(., fill=T, blank.lines.skip = T) %>%
    .[, c("id", "laterality", "visit", "measure") := tstrsplit(patientID, "_")]

## removeOutliers and restrict to healthy individuals
filt <- dt[id %in% healthyFinal] %>%
  .[! id %in% moorFieldsExclude[, eid]] %>%
  .[! patientID %in% outliers[,V1]] %>%
  .[, instances := paste(visit, measure, sep="_")]

sampsIDs <- filt[,.(id, patientID)]

fwrite(sampsIDs, file = "/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/data/scanIDs.txt", quote=F, col.names=F, sep = "\t")

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
    dt <- fread(file, col.names = c("patientID", "eye", "slice_index", "fovea_index", 0:255)) %>%
      .[, c("id", "laterality", "visit", "measure") := tstrsplit(patientID, "_")]

    fullScanMat <- dt[, c("slice_index", 0:255)] %>%
      as.matrix(., rownames="slice_index")


  # set -1 to NA
  fullScanMat[fullScanMat == -1] <- NA

  # summarise missingness for each scan
  naCountRow <- apply(fullScanMat, 1, function(x) sum(is.na(x)))
  naCountCol <- apply(fullScanMat, 2, function(x) sum(is.na(x)))

## calculate scan-wise summaries
  myMean <- mean(fullScanMat, na.rm = T)
  mySD <- sd(fullScanMat, na.rm = T)

  myMedian <- median(fullScanMat, na.rm = T)
  myMAD <- mad(fullScanMat, na.rm = T)


  # identify no. extreme values
  outliers3SDscan <- sum(!fullScanMat %between% c(myMean-3*mySD, myMean+3*mySD), na.rm=T)
  outliers4SDscan <- sum(!fullScanMat %between% c(myMean-4*mySD, myMean+4*mySD), na.rm=T)
  outliers5SDscan <- sum(!fullScanMat %between% c(myMean-5*mySD, myMean+5*mySD), na.rm=T)

  outliers3MADscan <- sum(!fullScanMat %between% c(myMedian-3*myMAD, myMedian+3*myMAD), na.rm=T)
  outliers4MADscan <- sum(!fullScanMat %between% c(myMedian-4*myMAD, myMedian+4*myMAD), na.rm=T)
  outliers5MADscan <- sum(!fullScanMat %between% c(myMedian-5*myMAD, myMedian+5*myMAD), na.rm=T)

  # data table to retun
    summ <- dt[, c("patientID", "eye", "id", "laterality", "visit", "measure")] %>%
      unique %>%
      .[, meanDepth := myMean] %>%
      .[, medianDepth := myMedian] %>%
      .[, minDepth := min(fullScanMat, na.rm = T)] %>%
      .[, maxDepth := max(fullScanMat, na.rm = T)] %>%
      .[, sdDepth := mySD] %>%
      .[, MADDepth := myMAD] %>%
      .[, totalNonMissing := sum(!is.na(fullScanMat))] %>%
      .[, totalMissing := sum(naCountRow)] %>%
      .[, meanMissingByRow := mean(naCountRow)] %>%
      .[, meanMissingByCol := mean(naCountCol)] %>%
      .[, outliers3SDscan := outliers3SDscan] %>%
      .[, outliers4SDscan := outliers4SDscan] %>%
      .[, outliers5SDscan := outliers5SDscan] %>%
      .[, outliers3MADscan := outliers3MADscan] %>%
      .[, outliers4MADscan := outliers4MADscan] %>%
      .[, outliers5MADscan := outliers5MADscan]

    }

    }) %>%
    rbindlist
  return(sampSumm)
}


## Get summaries for each scan.
cores <- 12
scansSumm <- mclapply(uniqueSamps, summlist, mc.cores = 12)

names(scansSumm) <- uniqueSamps
save(scansSumm, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/scansSummaries20221118.RData")


#load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/scansSummaries20220427.RData")



## Plot summaries
fields <- c("meanDepth", "medianDepth", "minDepth",  "maxDepth", "sdDepth", "MADDepth",
            "totalNonMissing", "totalMissing",  "meanMissingByRow", "meanMissingByCol",
            "outliers3MADscan", "outliers4MADscan", "outliers5MADscan")


lapply(fields, function(f) {

  plot <- scansSumm %>%
    rbindlist %>%
    select(.,  c("patientID", all_of(f)))

  ggplot(plot, aes_string(x = f)) +
    geom_histogram() +
    ggtitle(f) +
    scale_y_continuous(trans='log2')

    ggsave(file = paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/plots/summaries_",f,"_log.png"))

    ggplot(plot, aes_string(x = f)) +
      geom_histogram() +
      ggtitle(f)

    ggsave(file = paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/plots/summaries_",f,".png"))

} )


## Summarise available scans, by eye and instance
scansAvail <- mclapply(scansSumm , function(x) {

  if(nrow(x) > 0) {

  sc <- x[, c("id", "laterality", "visit", "measure") := tstrsplit(patientID, "_")]

   id <- sc[1,id]

   L_0 <-  sc[eye=="L" & visit==0] %>% nrow
   R_0 <- sc[eye=="R" & visit==0] %>% nrow
   L_1 <-  sc[eye=="L" & visit==1] %>% nrow
   R_1 <- sc[eye=="R" & visit==1] %>% nrow

  both_0 <- case_when(L_0 > 0 & R_0 > 0 ~ 1,
                     T ~ 0)
  both_1 <- case_when(L_1 > 0 & R_1 > 0 ~ 1,
                      T ~ 0)

  both_any <- case_when((L_0 + L_1) > 0 & (R_0 + R_1) > 0 ~ 1,
                        T ~ 0)

   out <- cbind(id, L_0, R_0, L_1, R_1, both_0, both_1, both_any) %>%
     as.data.table()

  return(out)

  }
  }, mc.cores = 12) %>%
  rbindlist


fwrite(scansAvail, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/scansAvailable20221118.txt", col.names=T)
