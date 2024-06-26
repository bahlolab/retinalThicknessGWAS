library(data.table)
library(magrittr)
library(dplyr)
library(here)

origResults <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/finalResultsEUR/gwSigLociSummary.csv") %>%
.[BonferroniSig == "Y"]

lociSNPs <-fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/finalResultsEUR/gwSigLociSNPs.csv") 
pixels <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/pixels.txt") 


## for each locus:
## extract the full list of SNPs for that locus
##  read in cojo results
## determine whether any locus SNPs are secondary signals
## check conditional analysis results for secondary signals

secondarySignals <- lapply(c(1:nrow(origResults)), function(locus) {

    print(paste0("Processing locus ", locus, " of ", nrow(origResults)))

    chr <- origResults[locus, CHR]
    pixel <- origResults[locus, pixel]
    SNPid <- origResults[locus, ID]

    locusSNPs <- lociSNPs[sentinelSNPID == SNPid] %>%
    .[!(locusSNPID == SNPid)] %>%
    .[, locusSNPID] 
     

    cojoFile <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/cojoAnalysis/output/chr", chr, "Pixel.", pixel, "_cojoOut.jma.cojo")
    condFile <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/cojoAnalysis/output/chr", chr, "Pixel", pixel, "_conditional_", SNPid, ".",pixel,".glm.linear")
    
    res <- fread(cojoFile) %>%
    .[SNP %in% locusSNPs] 

    if(nrow(res) > 0) {

        r2 <- lociSNPs[sentinelSNPID == SNPid  & locusSNPID %in% res[,SNP], .(locusSNPID, r2withSentinel)] %>%
        setnames(., "locusSNPID", "SNP") 

        condResult <- fread(condFile) %>%
          .[ID %in% res[,SNP]] %>%
          setnames(., "#POS", "POS") %>%
          .[, .(ID, POS, A1, BETA, SE, P)] %>%
          setnames(., c("ID", "POS_secondary_conditional", "A1", "BETA_secondary_conditional", "SE_secondary_conditional", "P_secondary_conditional")) %>%
          .[P_secondary_conditional < 5e-8 / 29041]


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
rsids <- fread("/stornext/Bioinf/data/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/rsIDsBonfSigSentinelsPixelwise.txt") %>%
    .[, CHR := ifelse(CHR == "X", "23", CHR) %>% as.integer]

newResults <- origResults %>%
    merge(., secondarySignals, by = c("ID")) %>%
    merge(., rsids, by = c("CHR", "POS")) 


fwrite(newResults, file ="/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/cojoAnalysis/output/bonfSigLociWithSecondary.csv", sep = ",")

newSentinels <- data.table(rsID = c(newResults[, rsID], newResults[, SNP_secondary]),
                            CHR = c(newResults[,CHR], newResults[, CHR]),
                            POS = c(newResults[, POS], newResults[, POS_secondary])) %>%
  na.omit %>%
  unique


## write vector as file
fwrite(newSentinels, file ="/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/cojoAnalysis/output/rsIDsBonfSigSentinelsSecondaryPixelwise.txt", sep = "\t")

# secondarySNPs <- newResults[!is.na(SNP_secondary), .(CHR, SNP_secondary)] %>%
#  unique

# secondaryResults <- lapply(c(1:nrow(secondarySNPs)), function(i) {
# snp <- secondarySNPs[i, SNP_secondary]
# chr <- secondarySNPs[i, CHR]

# print(paste("processing SNP:", snp, "on chr", chr))

# snpResults <- fread(here("secondarySignalResults", paste0(snp,".txt")))

# minP <- snpResults[P == min(P)]
# sigPix <-  snpResults[P < 5e-5] %>% nrow

# out <- minP %>%
#     .[, nPixelsSecondary := sigPix] %>%
#     .[, .(ID,  A1, BETA, SE, P, pixel, nPixelsSecondary)] %>%
#     setnames(., c("SNP_secondary", "A1_secondary", "BETA_topPixel_secondary", "SE_topPixel_secondary", "P_topPixel_secondary", "topPixel_Secondary", "nPixels_Secondary"))
# }) %>%
# rbindlist %>%
# .[newResults, on = "SNP_secondary"]
