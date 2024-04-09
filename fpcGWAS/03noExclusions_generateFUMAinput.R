library(data.table)
library(magrittr)
library(tidyverse)



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