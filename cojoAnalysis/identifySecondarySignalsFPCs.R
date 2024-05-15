library(data.table)
library(magrittr)
library(dplyr)
library(here)
library(stringr)
library(ieugwasr)



origResults <- fread("/stornext/Bioinf/data/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/sentinels/allChr_sentinel_clumpThresh0.001_withOverlap.txt") %>%
.[P < 5e-8/6]




## for each locus:
## extract the full list of SNPs for that locus
##  read in cojo results
## determine whether any locus SNPs are secondary signals
## check conditional analysis results for secondary signals

secondarySignals <- lapply(c(1:nrow(origResults)), function(locus) {

    print(paste0("Processing locus ", locus, " of ", nrow(origResults)))

    chr <- origResults[locus, CHR]
    fpc <- origResults[locus, FPC]
    SNPid <- origResults[locus, ID]

    locusSNPs <-  origResults[locus, SNPsInLocus] %>% 
        str_split(., ",") %>%  
        unlist 
     

    cojoFile <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/cojoAnalysis/output/chr", chr, ".fpc", fpc, "_cojoOut.jma.cojo")
    condFile <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/cojoAnalysis/output/chr", chr, "_conditional_", SNPid, ".fpc",fpc,".glm.linear")
    
    res <- fread(cojoFile) %>%
    .[SNP %in% locusSNPs] 

    if(nrow(res) > 0) {


    plinkPath <- "/vast/scratch/users/jackson.v/retThickness/GWAS/plink/plink2"
    bfilePath <-  paste0("/vast/scratch/users/jackson.v/retThickness/GWAS/geneticData/cleanedEURData/plinkBin/EUR_minMaf0.005_minInfo0.8_chr",chr)

     ld <- ld_matrix_local(c(SNPid, res[,SNP]), 
        bfile = bfilePath, 
        plink_bin = "/stornext/System/data/apps/plink/plink-1.9/bin/plink")

      ldSent <- ld[,colnames(ld) %like% SNPid]

      ldSNPs <- sapply(strsplit(names(ldSent), "_", fixed = TRUE),
       function(i) paste(head(i, -2), collapse = "_"))

      r2 <- data.table(SNP = ldSNPs, 
                              r2withSentinel = ldSent^2) %>%
                              .[SNP != SNPid]

        condResult <- fread(condFile) %>%
          .[ID %in% res[,SNP]] %>%
          setnames(., "#POS", "POS") %>%
          .[, .(ID, POS, A1, BETA, SE, P)] %>%
          setnames(., c("ID", "POS_secondary_conditional", "A1", "BETA_secondary_conditional", "SE_secondary_conditional", "P_secondary_conditional")) %>%
          .[P_secondary_conditional < 5e-8 / 6]


        if(nrow(condResult) > 0) {

                out <- res %>%
                .[, ID := SNPid] %>%
                .[r2, on = "SNP"] %>%
                .[condResult, on = c("SNP" = "ID", "refA" = "A1")] %>%
                .[, .(ID, SNP,   POS_secondary_conditional, refA, r2withSentinel, b, se, p,  BETA_secondary_conditional, SE_secondary_conditional, P_secondary_conditional)] %>%
                setnames(., c("ID", "SNP_secondary",  "POS_secondary", "A1_seconday", "r2_secondary", 
                    "BETA_secondary", "SE_secondary", "P_secondary",
                    "BETA_secondary_conditional", "SE_secondary_conditional", "P_secondary_conditional"))
        } else {    

        out <- data.table(ID = SNPid,
                            SNP_secondary = NA,
                            POS_secondary = NA,
                            A1_seconday = NA,
                            r2_secondary = NA,
                            BETA_secondary = NA,
                            SE_secondary = NA,
                            P_secondary = NA,
                            BETA_secondary_conditional = NA,
                            SE_secondary_conditional = NA,
                            P_secondary_conditional = NA)
        }                    

    } else {
        out <- data.table(ID = SNPid,
                            SNP_secondary = NA,
                            POS_secondary = NA,
                            A1_seconday = NA,
                            r2_secondary = NA,
                            BETA_secondary = NA,
                            SE_secondary = NA,
                            P_secondary = NA,
                            BETA_secondary_conditional = NA,
                            SE_secondary_conditional = NA,
                            P_secondary_conditional = NA)
    }

    return(out)
}) %>%
rbindlist




## merge origResults with secondarySignals, on "ID" and "A1".
rsids <- fread("/stornext/Bioinf/data/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/rsIDsBonfSigSentinelsFPC.txt") %>%
    .[, CHR := ifelse(CHR == "X", "23", CHR) %>% as.integer]

newResults <- origResults[, .(CHR, POS, ID,  A1, BETA, SE, P, FPC, clumpFPCs)] %>%
    merge(., secondarySignals, by = c("ID")) %>%
    merge(., rsids, by = c("CHR", "POS")) 


fwrite(newResults, file ="/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/cojoAnalysis/output/bonfSigFPCLociWithSecondary.csv", sep = ",")

newSentinels <- data.table(rsID = c(newResults[, rsID], newResults[, SNP_secondary]),
                            CHR = c(newResults[,CHR], newResults[, CHR]),
                            POS = c(newResults[, POS], newResults[, POS_secondary])) %>%
  na.omit %>%
  unique


## write vector as file
fwrite(newSentinels, file ="/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/cojoAnalysis/output/rsIDsBonfSigSentinelsSecondaryFPCs.txt", sep = "\t")

