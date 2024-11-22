library(data.table)
library(magrittr)
library(tidyverse)
library(ggmanh)
library(RColorBrewer)


chroms <- c(1:22, "X")
fpcs <- c(1:6)

####################
## rough code

## read in results and write out FUMA input file
lapply(fpcs, function(i) {
  
result <- lapply(chroms, function(chr) {
  
res <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/results/chr",chr,"/chr",chr,"EUR.fpc",i,".glm.linear") %>%
  fread %>%
  setnames(., "#CHROM", "CHR")

return(res)

}) %>%
rbindlist



FUMAResults <- result %>%
  .[, .(CHR, POS, REF, ALT, OBS_CT,  P)]
fwrite(FUMAResults, file =paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/FUMAresults/fpc",i,"_FUMA.txt"), sep="\t")

system(paste0("gzip -f /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/FUMAresults/fpc",i,"_FUMA.txt"),)
})


lapply(fpcs, function(i) {
  
hits <- lapply(chroms, function(chr) {

if(chr == "X") {
file <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/clumpedResults/chr",chr,"/chr",chr,"EUR.",i,".clumped") 

} else {
file <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/clumpedResults/chr",chr,"/chr",chr,"EUR.",i,"_thresh0.001_withOverlap.clumped") 
} 

if(file.exists(file)) {

res <- fread(file) %>%
  .[, .(SNP, CHR, BP)] %>%
  .[, CHR := ifelse(CHR==23, "X", CHR)]

return(res)

}

}) %>%
rbindlist



fwrite(hits, file =paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/FUMAresults/fpc",i,"_FUMAsentinels.txt"), sep="\t")

})




## read in results and write out for Yue
fpcResults <- lapply(c(1:8), function(i) {
  
result <- lapply(chroms, function(chr) {
  
res <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/results/chr",chr,"/chr",chr,"EUR.fpc",i,".glm.linear") %>%
  fread %>%
  setnames(., "#CHROM", "CHR")

return(res)

}) %>%
rbindlist %>%
  .[, .(CHR, POS, REF, ALT, A1_FREQ, BETA, SE, P)]
 return(result) 

}) %>%
rbindlist(., idcol = "fpc") %>%
dcast(., CHR + POS + REF + ALT + A1_FREQ ~ fpc, value.var = c("BETA", "SE", "P"))

fwrite(fpcResults, file ="/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/forYue/fpcs_genomeWide.txt", sep="\t")


lapply(fpcs, function(i) {
  
result <- lapply(chroms, function(chr) {
  
res <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/results/chr",chr,"/chr",chr,"EUR.fpc",i,".glm.linear") %>%
    fread %>%
    setnames(., "#CHROM", "CHR") %>%
    setnames(., "POS", "BP") %>%
    setnames(., "ID", "SNP") %>%
    setnames(., "REF", "A2") %>%
    setnames("A1_FREQ", "FREQ") %>%
    setnames(., "OBS_CT", "N") %>%
    return(res)
  
}) %>%
  rbindlist



outResults <- result %>%
  .[, .(CHR, BP, SNP, A1, A2, FREQ, BETA, SE, P, N)]
fwrite(outResults, file =paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/FUMAresults/fpc",i,"_CTG-VL.txt"), sep="\t")

})

## format for GWAS catalogue
## with these cols: .(chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, effect_allele_frequency, p_value, variant_id, n)

lapply(fpcs, function(i) {

print(i)  

result <- lapply(chroms, function(chr) {

res <- paste0("/stornext/Bioinf/data/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/results/chr",chr,"/chr",chr,"EUR.fpc",i,".glm.linear") %>%
  fread %>%
  setnames(., "#CHROM", "chromosome") %>%
  setnames(., "POS", "base_pair_location") %>%
  setnames(., "A1", "effect_allele") %>%
  setnames(., "REF", "other_allele") %>%
  setnames(., "A1_FREQ", "effect_allele_frequency") %>%
  setnames(., "BETA", "beta") %>%
  setnames(., "SE", "standard_error") %>%
  setnames(., "P", "p_value") 

  return(res)

}) %>%
rbindlist

## remove rows where alleles are not A, T, C, G
outResults <- result %>%
  .[chromosome == "X", chromosome := 23] %>%
 .[, .(chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, effect_allele_frequency, p_value)]

fwrite(outResults, file =paste0("/vast/scratch/users/jackson.v/retThickness/GWAS/FPCresultsGWAScat/FPC",i,"_GWAScatalogue.tsv"), sep="\t")
  
  })


## plot manhattan plot for each FPC using ggmanh
## save as pdf 
apply(fpcs, function(i) {

print(i)

result <- lapply(chroms, function(chr) {

print(chr)

res <- paste0("/stornext/Bioinf/data/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/results/chr",chr,"/chr",chr,"EUR.fpc",i,".glm.linear") %>%
  fread  %>%
  setnames(., "#CHROM", "CHR") %>%
  .[, .(CHR, POS, P)] 
return(res)
}) %>%
rbindlist

## set colours for manhattan plot
colours <- rep(brewer.pal(3,"PRGn")[c(1,3)], times = 12)


## set threshold for significant hits
sig <- 5e-8/6

## plot
plot <- manhattan_plot(result, 
          chr.colname = "CHR", 
          pos.colname = "POS",
          pval.colname = "P", 
          signif = sig, 
          chr.order = c(1:22, "X"),
          chr.col = colours[1:23],
          preserve.position = T,
          thin.n = 300,
          rescale = F)

ggsave(plot, 
       filename = paste0("/vast/scratch/users/jackson.v/retThickness/GWAS/FPCresultsManhattans/FPC",i,"_manhattan.pdf"),
       device = "pdf",
       height = 3, width = 8, units = "in")
})  


