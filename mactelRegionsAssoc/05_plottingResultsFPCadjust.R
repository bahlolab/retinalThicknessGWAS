#!/usr/bin/env Rscript

library(data.table)
library(magrittr)
library(tidyverse)
library(here)
library(corrplot)
library(gridExtra)

resultsDir <- "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/output/FPCadjust"

loci <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/rawData/macTelLoci.txt") %>%
    setnames(., c("Rsid", "chr", "BP", "from", "to"))
load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenoAssociations/rawData/Data_for_retinal_areas.RData")


scansDT  <-fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/rawData/scansUnadjustedFinal.csv", nrows=0)
pixels <- data.table(pixel = names(scansDT)[!names(scansDT) %in% c("patID", "visit", "eye", "sex", "age", "device", "meanRefErr")]) %>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert=TRUE)]


pixelWise <- lapply(loci[,Rsid], function(snp) {

snpResult <- lapply(pixels[,pixel], function(pix) {

   slice <- pixels[pixel==pix, y]

   print(paste(snp,pix))
  
   res <-  fread( paste0(resultsDir,"/",snp,"/",slice,"/",snp,"Pixel.",pix,".glm.linear")) %>%
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

  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/output/plots/",snp,"AssocsPixelwiseFPCadjust_",stat,".png"), width = 600, height = 600)
  print(plot)
  dev.off()

}

fwrite(snpResult, file = paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/output/",snp,"AssocsPixelwiseFPCadjust.txt"))

})



load("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/rawData/example_data.RData")

allLoci <- lapply(loci[,Rsid], function(snp) {

    result <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/output/",snp,"AssocsPixelwiseFPCadjust.txt") %>%
        fread %>%
        .[, .(pixel, y, x, BETA, P)] %>%
        .[, SNP := snp]

    # Plot the betas
    p1 <- ggplot(result)+
        geom_tile(aes(x=x,y=y,fill=BETA))+
        scale_fill_gradient2()+
    geom_path(aes(x=col,y=diag1),color="grey50")+
    geom_path(aes(x=col,y=diag2),color="grey50")+
    geom_path(aes(x=col,y=line1),color="grey50")+
    geom_path(aes(x=col,y=line2),color="grey50")+
    geom_path(aes(x=col,y=row,group=area),color="grey50",data = areas[areas$edtrs,],size = 0.5)+
    theme_bw()+theme(legend.position = "bottom")+ggtitle(paste(snp, "Beta"))+scale_y_reverse()

    # Plot the pvalues
    p2 <- ggplot(result)+
        geom_tile(aes(x=x,y=y,fill=P))+
        scale_fill_gradient2()+
    geom_path(aes(x=col,y=diag1),color="grey50")+
    geom_path(aes(x=col,y=diag2),color="grey50")+
    geom_path(aes(x=col,y=line1),color="grey50")+
    geom_path(aes(x=col,y=line2),color="grey50")+
    geom_path(aes(x=col,y=row,group=area),color="grey50",data = areas[areas$edtrs,],size = 0.5)+
    theme_bw()+theme(legend.position = "bottom")+ggtitle(paste(snp, "P-value"))+scale_y_reverse()


    plot <- gridExtra::grid.arrange(p1,p2,nrow=1)

    png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/output/plots/",snp,"AssocsPixelwiseFPCadjust_etdrs.png"), width = 1200, height = 600)
    print(plot)
    dev.off()

    return(result)
}) %>% 
rbindlist  


