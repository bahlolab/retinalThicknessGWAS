library(data.table)
library(magrittr)


## full sig results
dir <- "/vast/scratch/users/jackson.v/retThickness/fpcGWASnoExclusions/nonEURGWAS/results"

system("mkdir -p /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/AFRsentinels")
system("mkdir -p /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/CSAsentinels")

resultsAFR <- lapply(c(1:22, "X"), function(chr) {

print(paste(chr))


  FPCsResults <- lapply(1:6, function(i) {
    
    file <- paste0(dir,"/chr",chr,"/chr",chr,"_AFR.fpc",i,".glm.linear") 
    if(file.exists(file)) {
    fpcResult <-  fread(file) %>%
        setnames(., "#CHROM", "CHROM") %>%
        .[, .(CHROM, POS, ID, A1, A1_FREQ, BETA, SE, P)] %>%
        .[, FPC := i]

    } else {
        fpcResult <- NULL
    }

    return(fpcResult)

  }) %>%
  rbindlist

  return(FPCsResults)

}) %>%
rbindlist

fwrite(resultsAFR, file = paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/AFRsentinels/allChr_EURsentinels_AFRassocsAllFPCs.txt"))


resultsCSA <- lapply(c(1:22, "X"), function(chr) {

print(paste(chr))


  FPCsResults <- lapply(1:6, function(i) {
  
    file <- paste0(dir,"/chr",chr,"/chr",chr,"_CSA.fpc",i,".glm.linear") 
    if(file.exists(file)) {
    fpcResult <-  fread(file) %>%
        setnames(., "#CHROM", "CHROM") %>%
        .[, .(CHROM, POS, ID, A1, A1_FREQ, BETA, SE, P)] %>%
        .[, FPC := i]

    } else {
        fpcResult <- NULL
    }

    return(fpcResult)

  }) %>%
  rbindlist

  return(FPCsResults)

}) %>%
rbindlist

fwrite(resultsCSA, file = paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/CSAsentinels/allChr_EURsentinels_CSAassocsAllFPCs.txt"))

