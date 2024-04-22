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

secondarySignals <- lapply(c(1:nrow(origResults)), function(locus) {

    print(paste0("Processing locus ", locus, " of ", nrow(origResults)))

    chr <- origResults[locus, CHR]
    pixel <- origResults[locus, pixel]
    SNPid <- origResults[locus, ID]

    locusSNPs <- lociSNPs[sentinelSNPID == SNPid] %>%
    .[!(locusSNPID == SNPid)] %>%
    .[, locusSNPID] 
     

    file <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/cojoAnalysis/output/chr", chr, "Pixel.", pixel, "_cojoOut.jma.cojo")

    res <- fread(file) %>%
    .[SNP %in% locusSNPs] 

    if(nrow(res) > 0) {

        r2 <- lociSNPs[sentinelSNPID == SNPid  & locusSNPID %in% res[,SNP], .(locusSNPID, r2withSentinel)] %>%
        setnames(., "locusSNPID", "SNP") 

        out <- res %>%
        .[, ID := SNPid] %>%
        .[r2, on = "SNP"] %>%
        .[, .(ID, SNP,  refA, r2withSentinel, b, se, p)] %>%
        setnames(., c("ID", "SNP_secondary", "A1_seconday", "r2_secondary", "BETA_secondary", "SE_secondary", "P_secondary"))

    } else {
        out <- data.table(ID = SNPid,
                            SNP_secondary = NA,
                            A1_seconday = NA,
                            r2_secondary = NA,
                            BETA_secondary = NA,
                            SE_secondary = NA,
                            P_secondary = NA)
    }

    return(out)
}) %>%
rbindlist

## merge origResults with secondarySignals, on "ID" and "A1".

newResults <- origResults %>%
    merge(., secondarySignals, by = c("ID"))


secondarySNPs <- newResults[!is.na(SNP_secondary), .(SNP_secondary, CHR)] %>%
 unique

fwrite(secondarySNPs, file = here("secondarySignals.txt"))

# snp <- secondarySNPs[1, SNP_secondary]
# chr <- secondarySNPs[1, CHR]

# print(paste("processing SNP:", snp, "on chr", chr))

# snpResults <- lapply(c(1:nrow(pixels)), function(idx) {

#     if(idx %% 1000 == 0) {
#         print(paste0("Processing pixel ", idx, " of ", nrow(pixels)))
#     }

#     pix <-  pixels[idx, V1]
#     slice <- pixels[idx, V2]

#     file <- here("results", paste0("chr",chr), slice, paste0("chr",chr,"Pixel.",pix,".glm.linear.gz"))

#     ## fread, while grepping snp 
#     res <- fread(cmd=paste("zgrep", snp, file))
#     return(res)

# }) %>%
# rbindlist %>%
# setnames(. c("BP", "SNP", "A1", "b", "se", "Tstat", "p"))

# minP <- snpResults[p == min(p)]
# sigPix <-  snpResults[p < 5e-5] %>% nrow

# out <- minP %>%
#     .[, nPixelsSecondary := sigPix] %>%
#     .[, .(ID,  A1, b, se, p, nPixelsSecondary)] %>%
#     setnames(., c("SNP_secondary", "A1_secondary", "BETA_secondary", "SE_secondary", "P_secondary", "nPixels_Secondary"))
# )