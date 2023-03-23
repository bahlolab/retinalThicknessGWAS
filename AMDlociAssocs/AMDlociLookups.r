
library(data.table)
library(magrittr)
library(tidyverse)
library(patchwork)
library(stringr)
library(gridExtra)
library(grid)
library(gtable)
library(splyr)
library(RColorBrewer)


amd <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/AMDlociAssocs/rawData/Han_etal_2020_leadSNPs.csv") 

colnames(amd) <- gsub(" ", "_", colnames(amd)) 
colnames(amd) <- gsub("\\(|\\)", "", colnames(amd)) 


pixels <- fread("/vast/scratch/users/jackson.v/retThickness/GWAS/pixels.txt") %>%
setnames(., c("pixel", "y", "x"))

results <- lapply(c(1:22), function(chr) {

print(paste(chr))

dir <- paste0("/vast/scratch/users/jackson.v/retThickness/GWAS/GWsigResults/chr",chr)

chrResults <- lapply(c(1:119), function(slice) {


  file <- paste0(dir,"/chr",chr,"Slice",slice,"_5e-5Sig.txt")

  sliceResults <- fread(file) %>%
    setnames(., "#POS", "POS") %>%
    .[ID %in% amd[,SNP]]
  return(sliceResults)

}) %>%
rbindlist

return(chrResults)
}) %>%
rbindlist



sentinelIDs <- results[,ID] %>% unique

lapply(sentinelIDs, function(snp) {

snpOut <- ifelse(snp %like% ":", str_replace(snp, ":", "_"), snp)

snpResult <- results[ID == snp] %>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert = T)] %>%
  .[!is.na(y)] %>%
  .[, log10P := log(P, 10)]


statsPlots <- lapply(c("BETA", "T_STAT", "log10P"), function(stat) {
    plot <- ggplot(snpResult) +
    geom_tile(aes_string(x = "x", y = "y", fill = stat)) +
    # geom_path(aes(x=col,y=row,color=area),data = areas, size = 0.5) +
    scale_fill_gradient2() +
    scale_y_reverse() +
    theme_bw() +
    theme(legend.position = "bottom")+
    ggtitle(paste(snp,stat, sep = " - "))

  return(plot)
})

  png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/AMDlociAssocs/output/plots/assocsPixelwise_",snpOut,".png"), width = 1800, height = 600)
  reduce(statsPlots , `+`) %>%
  print
  dev.off()

})



fpcResults <-  lapply(c(1:22), function(chr) {

print(paste(chr))

chrResults <- lapply(c(1:8), function(fpc) {
  
  # print(fpc)
  # slice <- 64
  file <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/results/chr",chr,"/chr",chr,"EUR.fpc",fpc,".glm.linear")
  
  fpc <- fread(file) %>%
    setnames(., "#CHROM", "CHR") %>%
    .[ID %in% amd[,SNP]] %>%
    .[, FPC := fpc]
  
  return(fpc)
  
}) %>%
  rbindlist

return(chrResults)
}) %>%
rbindlist



fpcLogP <- fpcResults %>%
  .[, log10P := (-1)*log(P, 10)] %>%
  .[, sig := case_when(P <= 5E-5 & P > 5E-8 ~ "<5E-5",
                        P <= 5E-8 ~ "<5E-8",
                        T ~ "")]

order <- fpcLogP[order(POS, decreasing = T),ID] %>%
unique

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/AMDlociAssocs/output/plots/AMDsentinels_fpcResults.png"), width = 1200, height = 1200)
ggplot(fpcLogP, aes(y = ID, x = FPC)) +
    geom_tile(aes(fill = log10P)) +
    geom_text(aes(y=ID,  x = FPC, label=sig)) +
    scale_fill_gradient2() +
    # scale_y_reverse()  + 
    scale_y_discrete(limits = order) +
    scale_x_continuous(breaks=seq(1,10,1)) +
    theme_bw() +
    theme(legend.position = "bottom")
dev.off()




minP <-  results[, minP := min(.SD), .SDcols = "P", by = "ID"] %>%
    .[,.(ID, minP)] %>%
    unique 
  
 minPamd <-   amd[, .(SNP, Nearest_Gene, P_MTAG)] %>%
    setnames(., c("ID", "gene", "P_AMD")) %>%
  minP[., on = "ID"] %>%
  .[, retThickness := (-1)*log(minP, 10)] %>%
  .[, AMD := (-1)*log(as.numeric(P_AMD), 10)] %>%
  .[, gwSig := ifelse(minP < 5E-8, 1, 0)] %>%
  .[, .(ID, gene, retThickness, AMD, gwSig)] %>%
  melt(., id.vars = c("ID", "gene", "gwSig"), variable.name = "Analysis", value.name = "log10P") %>% 
  .[, log10P := if_else(Analysis == "retThickness", -log10P, log10P)] %>%
  .[, col := case_when(gwSig==1 ~ 2,
      gwSig==0 ~ 1,
      T ~ 0)]

delete_col <- function(x, pattern) {
  t <- x$layout %>% 
    filter(str_detect(name, pattern)) %>% 
    pull(l)

  x <- gtable_filter(x, pattern, invert = TRUE)

  x$widths[t] <- unit(0, "cm")

  x
}

order <- minPamd[Analysis == "AMD"] %>%
.[order(log10P)] %>%
.[,ID]

labs <- minPamd[Analysis == "AMD"] %>%
.[order(log10P)] %>%
.[,gene]

## plot 1


p1 <- minPamd %>% 
  ggplot(aes(ID, log10P)) + 
  facet_wrap(~ Analysis, scales = "free_x") + 
  geom_col() +
  coord_flip() +
  theme(axis.text.y = element_text(hjust = 0.5, margin = margin(0, 0, 0, 0)),
        axis.ticks.length = unit(0, "pt")) +
  labs(x = NULL)  +
  scale_x_discrete(limits = order, labels = labs) + 
  theme_minimal() + 
   scale_fill_manual(values =  c("#666666"))

p1_g <- ggplotGrob(p1)

p1_g$widths[7] <- p1_g$widths[4] + unit(0.5, "cm")

p1g_axis <- gtable_filter(p1_g, "axis-l-1-1") 

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/AMDlociAssocs/output/plots/AMD_retThickness_compare.png"), width = 1200, height = 1200)
p1_g %>% 
  gtable_add_grob(p1g_axis, l = 7, t = 8, name = "Locus") %>% # add the axis to the middle
  delete_col("axis-l-1-1") %>% # delete the original axis
  gtable_add_grob(textGrob("Locus", gp = gpar(fontsize = 11)), l = 7, t = 7) %>% # add the top label
  grid.draw() # draw the result
dev.off()



p2 <- minPamd %>% 
  ggplot(aes(ID, log10P, fill = as.factor(col) )) + 
  facet_wrap(~ Analysis, scales = "free_x") + 
  geom_col() +
  coord_flip() +
  theme(axis.text.y = element_text(hjust = 0.5, margin = margin(0, 0, 0, 0)),
        axis.ticks.length = unit(0, "pt")) +
  labs(x = NULL)  +
  scale_x_discrete(limits = order, labels = labs) + 
  theme_minimal() +
   guides(fill = "none") + 
   scale_fill_manual(values =  c( "#666666", "#D95F02", "#7570B3"))

p2_g <- ggplotGrob(p2)

p2_g$widths[7] <- p2_g$widths[4] + unit(0.5, "cm")

p2g_axis <- gtable_filter(p2_g, "axis-l-1-1") 

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/AMDlociAssocs/output/plots/AMD_retThickness_compare_retThickSig.png"), width = 1200, height = 1200)
p2_g %>% 
  gtable_add_grob(p2g_axis, l = 7, t = 8, name = "Locus") %>% # add the axis to the middle
  delete_col("axis-l-1-1") %>% # delete the original axis
  gtable_add_grob(textGrob("Locus", gp = gpar(fontsize = 11)), l = 7, t = 7) %>% # add the top label
  grid.draw() # draw the result
dev.off()


# plot 3
order <- minPamd[Analysis == "AMD" & !is.na(gwSig)] %>%
.[order(log10P)] %>%
.[,ID]

labs <- minPamd[Analysis == "AMD" & !is.na(gwSig)] %>%
.[order(log10P) ] %>%
.[,gene]

p3 <- minPamd[!is.na(gwSig)] %>% 
  ggplot(aes(ID, log10P, fill = as.factor(col) )) + 
  facet_wrap(~ Analysis, scales = "free_x") + 
  geom_col() +
  coord_flip() +
  theme(axis.text.y = element_text(hjust = 0.5, margin = margin(0, 0, 0, 0)),
        axis.ticks.length = unit(0, "pt")) +
  labs(x = NULL)  +
  scale_x_discrete(limits = order, labels = labs) + 
  theme_minimal() +
   guides(fill = "none") + 
   scale_fill_manual(values =  c("#D95F02", "#7570B3"))

p3_g <- ggplotGrob(p3)

p3_g$widths[7] <- p3_g$widths[4] + unit(0.5, "cm")

p3g_axis <- gtable_filter(p3_g, "axis-l-1-1") 

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/AMDlociAssocs/output/plots/AMD_retThickness_compare_retThickSigOnly.png"), width = 1200, height = 600)
p3_g %>% 
  gtable_add_grob(p3g_axis, l = 7, t = 8, name = "Locus") %>% # add the axis to the middle
  delete_col("axis-l-1-1") %>% # delete the original axis
  gtable_add_grob(textGrob("Locus", gp = gpar(fontsize = 11)), l = 7, t = 7) %>% # add the top label
  grid.draw() # draw the result
dev.off()
