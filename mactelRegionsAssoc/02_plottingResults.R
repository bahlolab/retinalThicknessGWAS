#!/usr/bin/env Rscript

library(data.table)
library(magrittr)
library(tidyverse)
library(here)
library(corrplot)
library(gridExtra)
library(funData)
library(patchwork)
library(purrr)

resultsDir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/output"

loci <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/rawData/macTelLoci.txt") %>%
    setnames(., c("Rsid", "chr", "BP", "from", "to"))
load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenoAssociations/rawData/Data_for_retinal_areas.RData")


scansDT  <-fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/rawData/scansUnadjustedFinal.csv", nrows=0)
pixels <- data.table(pixel = names(scansDT)[!names(scansDT) %in% c("patID", "visit", "eye", "sex", "age", "device", "meanRefErr")]) %>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert=TRUE)]

# finished <- c("rs662602", "rs17421627", "rs146953046")

pixelWise <- lapply(loci[,Rsid], function(snp) {
# pixelWise <- lapply(finished, function(snp) {
# pixelWise <- lapply(loci[!(Rsid %in% finished),Rsid], function(snp) {


pos <- loci[Rsid==snp, BP]

snpID <- paste0(resultsDir,"/",snp,"/1/",snp,"Pixel.1_39.glm.linear") %>%
    fread %>%
    .[POS==pos] %>%
    .[,ID]

colNames <- paste0(resultsDir,"/",snp,"/1/",snp,"Pixel.1_39.glm.linear") %>%
    fread(., nrows=0)

snpResult <- lapply(pixels[,pixel], function(pix) {

   slice <- pixels[pixel==pix, y]

   print(paste(snp,pix))
  
   res <-  fread(cmd = paste0("grep ",snpID," ", resultsDir,"/",snp,"/",slice,"/",snp,"Pixel.",pix,".glm.linear")) %>%
    setnames(., names(colNames)) %>%
    .[, pixel := pix]
  
   return(res) 
  }) %>%
  rbindlist() %>%
  .[, log10P := (-1)*log10(P)] %>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert=TRUE)]


for(stat in c("BETA", "T_STAT", "log10P")) {

  plot <- ggplot(snpResult) +
    geom_tile(aes_string(x = "x", y = "y", fill = stat)) +
    scale_fill_gradient2() +
    scale_y_reverse() +
    theme_bw() +
    theme(legend.position = "bottom")+
    ggtitle(paste(snp,stat, sep = " - "))




  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/output/plots/",snp,"AssocsPixelwise_",stat,".png"), width = 600, height = 600)
  print(plot)
  dev.off()

}

fwrite(snpResult, file = paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/output/",snp,"AssocsPixelwise.txt"))

})

fpcs <- c(1:50)


load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/rawData/example_data.RData")


## update rsis for PHGDH rare SNP (unsure why wrong...?)
loci <- loci[BP==120265444, Rsid := "rs532303"]

allLoci <- lapply(loci[,Rsid], function(snp) {

    result <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/output/",snp,"AssocsPixelwise.txt") %>%
        fread %>%
        .[, .(pixel, y, x, BETA, P)] %>%
        .[, SNP := snp]

    # Plot the betas
    p1 <- ggplot(result)+
        geom_tile(aes(x=x,y=y,fill=BETA))+
        scale_fill_gradient2()+
      geom_path(aes(x=col,y=diag1),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=diag2),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=line1),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=line2),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=row,group=area),color="grey50",data = areas[areas$edtrs,],size = 0.5)+
    theme_bw()+theme(legend.position = "bottom")+ggtitle(paste(snp, "Beta"))+scale_y_reverse()

    # Plot the pvalues
    p2 <- ggplot(result)+
        geom_tile(aes(x=x,y=y,fill=P))+
        scale_fill_gradient2()+
      geom_path(aes(x=col,y=diag1),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=diag2),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=line1),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=line2),color="grey50", data=res[res$variable=="logFC",])+
      geom_path(aes(x=col,y=row,group=area),color="grey50",data = areas[areas$edtrs,],size = 0.5)+
    theme_bw()+theme(legend.position = "bottom")+ggtitle(paste(snp, "P-value"))+scale_y_reverse()



    png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/output/plots/",snp,"AssocsPixelwise_etdrs.png"), width = 1200, height = 600)
    gridExtra::grid.arrange(p1,p2,nrow=1)
    dev.off()

    return(result)
}) %>% 
rbindlist  

 mat <- allLoci %>% 
    dcast(., pixel ~ SNP, value.var = "BETA") %>%
    as.matrix(., rownames = "pixel") %>%
     # abs() %>%
    na.omit 
  
corrMat <- cor(mat)

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/output/plots/macTelLoci_corrHeatmap.png")
corrplot(corrMat)
dev.off()


## fpc results
## read in results 
fpcs <- c(1:50)

fpcResults <- lapply(loci[,Rsid], function(snp) {

    chr <- loci[Rsid == snp,  chr]

    snpResult <- lapply(fpcs, function(i) {
    
        res <- paste0("/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/fpcGWAS/results/chr",chr,"/chr",chr,".fpc",i,".glm.linear") %>%
        fread %>%
        setnames(., "#CHROM", "CHR") %>%
        .[ID == snp]
        return(res)

    }) %>%
    rbindlist(idcol = "fpc") 
}) %>%
rbindlist() %>%
.[, log10P := (-1)*log10(P)] %>%
.[, gwSig := ifelse(P < 5e-8, 1, 0)]


for(stat in c("BETA", "T_STAT", "log10P", "gwSig")) {

  plot <- ggplot(fpcResults[fpc <= 25]) +
    geom_tile(aes_string(x = "ID", y = "fpc", fill = stat)) +
    scale_fill_gradient2() +
    scale_y_reverse() +
    theme_bw() +
    theme(legend.position = "bottom")+
    ggtitle(paste("MacTel Loci - fpcGWAS results",stat, sep = " - "))



  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/output/plots/MacTelLoci_fpcGWASResults_",stat,".png"), width = 600, height = 600)
  print(plot)
  dev.off()

}


plot <- ggplot(fpcResults[fpc <= 25]) +
  geom_tile(aes(x = ID, y = fpc, fill = T_STAT)) +
  geom_text(aes(x = ID, y = fpc, label = signif(P, 3)), data = fpcResults[fpc <= 25 & P < 5e-5]) +
  scale_fill_gradient2() +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "bottom")+
  ggtitle("MacTel Loci - fpcGWAS results")

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/output/plots/MacTelLoci_fpcGWASResults.png"), width = 600, height = 600)
print(plot)
dev.off()




load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenoExploratory/working/cleanedScansFPCA_20220927.RData")
scans <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/processedData/scansUnadjustedFinal.csv")
pixels <-  names(scans)[!names(scans) %in% c("patID", "eye", "visit", "sex", "age", "device", "meanRefErr")]



fpcPlots <- lapply(c(1:20), function(i) {

   fpc <- pca$functions[[1]] %>%
    funData::as.data.frame(.) %>%
    as.data.table %>%
    .[obs==i]

  plot <- ggplot(fpc) +
    geom_tile(aes(y = argvals1, x = argvals2, fill = X)) +
    scale_fill_gradient2() +
    scale_y_reverse() +
    theme_bw() +
    theme(legend.position = "none")
  # theme(legend.position = "bottom")

  return(plot)

})

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/generateFPCs/output/fpcFunctions.png", width = 1500, height = 1200)
reduce(fpcPlots , `+`) %>%
  print
dev.off()

