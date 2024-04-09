# BiocManager::install("ggmanh")
library(data.table)
library(magrittr)
library(ggmanh)

dataDir <- "/vast/scratch/users/jackson.v/retThickness/GWAS"
outDir <- "/stornext/Bioinf/data/lab_bahlo/projects/misc/retinalThickness/GWAS/plots/manhattans"

## set i as slice number
i <- 1

chroms <- c(1:22, "X")

## results for slice
sliceResult <- lapply(chroms, function(chr) {

print(paste("Reading results for chromosome", chr))    

res <- paste0(dataDir,"/chr",chr,"/slice",i,"_result.txt.gz") %>%
fread()

return(res)

}) 
names(sliceResult) <- chroms

## get pixels for slice
pixels <- grep(paste0("^",i,"_.*_P$"), names(sliceResult[[1]]), value = TRUE) %>%
    sub(paste0("^",i,"_(\\d+)_P$"), "\\1", .)

lapply(pixels, function(pix) {

## extract genome-wide result for pixel
cols <- c("POS", paste(i, pix, "P", sep = "_"))

pixResult <- lapply(chroms, function(chr) {

    sliceResult[[chr]][, ..cols] %>%
    setnames(., c("POS", "P")) %>%
    .[, CHR := chr] %>%
    return

}) %>%
rbindlist %>%
.[, CHR := as.factor(CHR)]

print(paste("Plotting pixel", pix))    

## manhattan plot
pixPlot <- manhattan_plot(x = pixResult, 
                        pval.colname = "P", 
                        chr.colname = "CHR", 
                        pos.colname = "POS", 
                        chr.order = chroms,
                        plot.title = paste0("Pixel ",i,"_",pix),
                        thin.n = 500,
                        signif = c(1.72e-12, 5e-08))
                        
png(paste0(outDir,"/",i,"/Pixel",i,"_",pix,"_manhattan.png"), width = 1000, height = 500)
print(pixPlot)
dev.off()

rm(pixResult)
rm(pixPlot)

})
