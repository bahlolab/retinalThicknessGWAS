library(data.table)
library(magrittr)
library(tidyverse)

chr <- 22
dir <- paste0("/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/results/chr",chr)
sentinels <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/chr",chr,"sentinels.csv") %>%
 fread %>%
 setnames(., "#POS", "POS")

 lapply(c(1:nrow(sentinels)), function(i) {

    SNP <- sentinels[i,ID]  

    pix <- sentinels[i, pixel]
    slice <- str_split(pix, "_") %>%
        unlist %>% 
        .[1]

    pos <- sentinels[i, POS]
    posMin <- pos - 1000000
    posMax <- pos + 1000000

    results <- paste0(dir,"/",slice,"/chr",chr,"Pixel.",pix,".glm.linear.gz") %>% 
    fread %>%
    setnames(., "#POS", "POS") %>%
    .[POS %between% c(posMin, posMax)] %>%
    .[, .(ID, P)]

    # locuszoom.commands <- data.table(snp=SNP,
    #                              chr=NA,
    #                              start=NA,
    #                              end=NA,
    #                              flank="500kb",
    #                              run="yes",
    #                              arguments="showAnnot=T showRecomb=T")


    fwrite(results, file = paste0("/vast/scratch/users/jackson.v/retThickness/GWAS/regionPlots/chr",chr,"/",SNP,"_METAL.txt"), sep = "\t")
    # fwrite(locuszoom.commands, file = "", sep = "\t")
 })

