library(data.table)
library(magrittr)

chr <- 1 
slice <- 1
pix <- "1_64"

# res <- lapply(1:22, function(chr) {
#   paste0("/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/results/chr",chr,"/",slice,"/chr",chr,"Pixel.",pix,".glm.linear.gz") %>%
#   fread(., select = c("ID", "P")) %>%
#   .[P < 5E-5] 
# })  %>%
#   rbindlist

pixels <- fread("/vast/scratch/users/jackson.v/retThickness/GWAS/pixels.txt") %>%
  setnames(., c("pixel", "y", "x"))  

nPix <- nrow(pixels)

nSNPs <- lapply(1:22, function(chr) {
  paste0("/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/results/chr",chr,"/",slice,"/chr",chr,"Pixel.",pix,".glm.linear.gz") %>%
    fread(., select = c("ID"), nrows = 5) %>%
    nrow 
})  %>%
  Reduce(sum, .)

n <- as.numeric(nPix)*as.numeric(nSNPs)


## full sig results
res <- lapply(1:22, function(chr) {
  print(paste("chr",chr))
dir <- paste0("/vast/scratch/users/jackson.v/retThickness/GWAS/GWsigResults/chr",chr)

results <- lapply(c(1:119), function(slice) {
  
  print(paste(slice))
  
  file <- paste0(dir,"/chr",chr,"Slice",slice,"_5e-5Sig.txt")
  
  sliceResults <- fread(file, select = c("ID", "P")) %>%
    .[P < 5E-5] 
  
  return(sliceResults)
  
}) %>%
  rbindlist
}) %>%
  rbindlist


p <- res[,P]
nm <- res[,ID]
p <- as.numeric(p)
p0 <- setNames(p, nm)

nna <- !is.na(p)
p <- p[nna]
lp <- length(p)

i <- lp:1L
o <- order(p, decreasing = TRUE)
ro <- order(o)


splits <- seq(1L, n, by = 1e6)

q <- lapply(1:length(splits), function(x) {
  
  start <- splits[x]
  if(x < length(splits)) {
    
    end <- splits[x+1]-1
  } else {
    end <- n
  }
  
  sum(1/(start:end)) %>% return()
      
}) %>%
  unlist %>%
  sum

p0 <- pmin(1, cummin(q * n/i * p[o]))[ro]

resultsOut <- res[, pAdj := p0]

fwrite(resultsOut, file = "/vast/scratch/users/jackson.v/retThickness/GWAS/GWsigResults/adjustedPvals.txt")
