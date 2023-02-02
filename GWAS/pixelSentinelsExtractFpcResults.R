

fpcResults <- lapply(c(1:10), function(fpc) {
  
  # print(fpc)
  # slice <- 64
  file <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/results/chr",chr,"/chr",chr,"EUR.fpc",fpc,".glm.linear")
  
  fpcResults <- fread(file) %>%
    setnames(., "#CHROM", "CHR") %>%
    .[, FPC := fpc]
  
  return(fpcResults)
  
}) %>%
  rbindlist


pixelSentinels <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/chr",chr,"sentinels_clumpThresh0.001.txt") %>%
  fread(.)

fpcLogP <- fpcResults[ID %in% pixelSentinels[,ID]] %>%
  .[, log10P := (-1)*log(P, 10)] %>%
  .[, sig := case_when(P < 5E-5 ~ "<5E-5",
                        P < 5E-8 ~ "<5E-8",
                        T ~ "")]

order <- fpcLogP[order(POS, decreasing = T),ID] %>%
unique

png(paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/plots/chr",chr,"/PixelwiseSentinels_fpcResults.png"), width = 1200, height = 1200)
ggplot(fpcLogP, aes(y = ID, x = FPC)) +
    geom_tile(aes(fill = log10P)) +
    geom_text(aes(y=ID,  x = FPC, label=sig)) +
    scale_fill_gradient2() +
    # scale_y_reverse()  + 
    scale_y_discrete(limits = order) +
    scale_x_continuous(breaks=seq(1,10,1)) +
    theme_bw() +
    theme(legend.position = "bottom")+
    ggtitle(paste("Chromosome",chr))
dev.off()
