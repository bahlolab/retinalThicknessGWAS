#!/usr/bin/env Rscript


library(data.table)
library(magrittr)
library(tidyverse)
library(RColorBrewer)
library(stringr) 
library(pheatmap)

corrs <- lapply(c(1:6), function(i) {
  paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/misc/genCorrs/fpc",i,"_ldsc.csv") %>%
    fread
})

sigTraits <- lapply(corrs, function(dt) {
  
  dt[as.numeric(p) < 0.01, p2] 
  
}) %>%
  unlist 

sigCorrs <- lapply(corrs, function(dt) {
  
  dt[p2 %in% sigTraits] 
  
  }) %>%
  rbindlist(., idcol = "FPC") %>%
  .[, sig := case_when(as.numeric(p) < 0.0005 ~ "***",
                       as.numeric(p) < 0.005 ~ "**",
                       as.numeric(p) < 0.05 ~ "*",
                       T ~ "")] %>%
  .[order(as.numeric(p))] %>%
  .[, rg :=  as.numeric(rg)] 




corrs <- lapply(c(1:6), function(i) {
  paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/misc/genCorrs/fpc",i,"_ldsc.csv") %>%
    fread
}) %>%
  rbindlist(., idcol = "FPC") %>%
  .[!str_detect(phenotype, "work|Work|job|Job|pollution")]
  

sigTraits <- corrs %>%
  .[order(as.numeric(p))] %>%
  .[1:50, p2]

sigCorrs <- corrs[p2 %in% sigTraits] %>%
  .[, sig := case_when(as.numeric(p) < 0.0005 ~ "***",
                       as.numeric(p) < 0.005 ~ "**",
                       as.numeric(p) < 0.05 ~ "*",
                       T ~ "")] %>%
  .[order(as.numeric(p))] %>%
  .[, rg :=  as.numeric(rg)] 

colours <- rev(brewer.pal(9,"RdBu"))

clust <- sigCorrs %>%
  dcast(., phenotype ~ FPC , value.var = "rg") %>%
  as.matrix(., rownames = "phenotype") %>%
  dist() %>%
  hclust()

sigCorrs %>% 
  ggplot(., aes(x = as.factor(FPC), y = phenotype)) +
  geom_raster(aes(fill = as.numeric(rg)), color = "white", lwd = 0.3,linetype = 1) +
  geom_text(aes(label = sig))  +
  labs(y = "Traits", title = "", x = "fPC") +
  theme(legend.position = "top", legend.key.width= unit(0.8, 'cm')) +
  scale_fill_gradientn(colors=colours, limits = c(-0.6, 0.6), name = "Genetic Correlation (rg)" ) +
  scale_y_discrete(limits = clust$labels[clust$order])
  
  
resultsDT <- corrs %>%
dcast(., phenotype + Cohort + URL + Notes ~ FPC , value.var = c("rg", "p"))
      
 fwrite(resultsDT, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/genCorrelations.csv")     
 