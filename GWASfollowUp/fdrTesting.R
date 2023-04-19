library(data.table)
library(magrittr)


## full sig results
res <- lapply(1:22, function(chr) {
  print(paste("chr",chr))
  fpcResults <- lapply(c(1:8), function(fpc) {
    
    file <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/results/chr",chr,"/chr",chr,"EUR.fpc",fpc,".glm.linear")
    
    fpcResults <- fread(file) %>%
      .[,.(ID, P)] %>%
      .[, FPC := fpc]
    
    return(fpcResults)
    
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

# fwrite(resultsOut, file = "/vast/scratch/users/jackson.v/retThickness/GWAS/GWsigResults/adjustedPvals.txt")

test <- resultsOut %>%
  .[, padjBYtest := p.adjust(res[,P], method = "BY")] %>%
  .[, padjBHtest := p.adjust(res[,P], method = "BH")] %>%
  .[, padjBontest := p.adjust(res[,P], method = "bonferroni")]



testRes <- test[P < 5E-5] %>%
  .[, padjBYtestAlt := p.adjust(P, method = "BY", n=nrow(res))] %>%
  .[, padjBHtestAlt := p.adjust(P, method = "BH", n=nrow(res))] %>%
  .[, padjBontestAlt := p.adjust(P, method = "bonferroni", n=nrow(res))]







# chr <- 1 
# slice <- 1
# pix <- "1_64"

# res <- lapply(1:22, function(chr) {
#   paste0("/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/results/chr",chr,"/",slice,"/chr",chr,"Pixel.",pix,".glm.linear.gz") %>%
#   fread(., select = c("ID", "P")) %>%
#   .[P < 5E-5] 
# })  %>%
#   rbindlist

pixels <- fread("/vast/scratch/users/jackson.v/retThickness/GWAS/pixels.txt") %>%
  setnames(., c("pixel", "y", "x"))  

nPix <- nrow(pixels)
slice <- 1
pix <- "1_64"

nSNPs <- lapply(c(1:22, "X"), function(chr) {
  paste0("/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/results/chr",chr,"/",slice,"/chr",chr,"Pixel.",pix,".glm.linear.gz") %>%
    fread(., select = c("ID")) %>%
    nrow 
})  %>%
  Reduce(sum, .)

n <- as.numeric(nPix)*as.numeric(nSNPs)


## full sig results
# fullRes <- lapply(1:22, function(chr) {
#   print(paste("chr",chr))
# dir <- paste0("/vast/scratch/users/jackson.v/retThickness/GWAS/GWsigResults/chr",chr)
# 
# sliceResults <- lapply(c(1:119), function(slice) {
#   
#   print(paste(slice))
#   
#   pixResults <- lapply(pixels[y==slice, pixel] , function(pix) {
#     
#   file <-  paste0("/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/results/chr",chr,"/",slice,"/chr",chr,"Pixel.",pix,".glm.linear.gz")
#   
#   res <- fread(file, select = c("ID", "P")) %>%
#     .[P < 5E-5] %>%
#     .[, pix := pix]
#   return(res)
#   }) %>%
#   rbindlist
#   
#   return(pixResults)
#   
# }) %>%
#   rbindlist
# 
# return(sliceResults)
# }) %>%
#   rbindlist

fullRes <- lapply(c("X", 22:1), function(chr) {
  print(paste("chr",chr))

  chrResult <- paste0("/vast/scratch/users/jackson.v/retThickness/GWAS/FDR/chr",chr,".txt") %>%
  fread(., header=F) %>%
    setnames(., c("ID", "P", "pix"))
  return(chrResult)
}) %>%
  rbindlist

p <- fullRes[,P]
nm <- fullRes[,ID]
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

fullRes <- fullRes[, pAdj := p0]

fullRes <- fullRes %>%
  # .[, padjBYtest := p.adjust(P, method = "BY", n=n)] %>%
  .[, padjBHtest := p.adjust(P, method = "BH", n=n)] %>%
  .[, padjBontest := p.adjust(P, method = "bonferroni", n=n)]

fwrite(fullRes, file = "/vast/scratch/users/jackson.v/retThickness/GWAS/FDR/pixelWiseAdjustedPvals.txt")
