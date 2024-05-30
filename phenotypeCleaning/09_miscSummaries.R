library(data.table)
library(magrittr)
library(dplyr)
library(rlist)

## raw scan data
scansDT <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/scansWideFormat.csv")

pixels <-  names(scansDT)[!names(scansDT) %in% c("patID", "eye", "visit", "scan",   "refErr")]

## missing data points for pixels only
scansMissing <- is.na(scansDT[, ..pixels]) %>%
    cbind(scansDT[, c("patID", "eye", "visit", "scan",   "refErr")], .) 

rm(scansDT)

## read in imputed scans
bestScansDT <- fread("/stornext/Bioinf/data/lab_bahlo/projects/misc/retinalThickness/phenoExploratory/bestScansIDT20221118.txt")

imputedScans <- lapply(c(1:46), function(chunk) {
  paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/scansImputedWideFormat_chunk",chunk,".csv") %>%
  fread %>%
    .[scan %in% bestScansDT[,scan]]
  }) %>% rbindlist 

rm(bestScansDT)

## for each pixel, obtain mean from corresponding column in imputedScans data.table, stratified by whether corresponding pixel colum in scansMissing is TRUE or FALSE

## In imputedScans, set values to missing, if corresponding value in scansMissing is TRUE
setkey(imputedScans, patID, eye, visit, scan, refErr)
setkey(scansMissing, patID, eye, visit, scan, refErr)

# create new DT missingvalues, with clumns pixels, where value is from imputedScans if scansMissing is TRUE, otherwise NA
missingnessMeans <- lapply(pixels, function(p) {

    # print every 100th pixel
    if(which(pixels==p) %% 100 == 0) {
        print(paste("processing pixel", which(pixels==p) ,"/", length(pixels)))
    }

    missing <- scansMissing[, ..p] %>% unlist

    missingMean <- imputedScans[missing, ..p] %>%
        unlist %>%
        mean(na.rm = TRUE)

    nonMissingMean <- imputedScans[!missing, ..p] %>%
        unlist %>%
        mean(na.rm = TRUE)

    ## test for difference in means between missing and non-missing
    pDiff <- t.test(imputedScans[missing, ..p] %>% unlist, imputedScans[!missing, ..p] %>% unlist) %$% p.value    

    out <- data.table(pixel = p,
                     missingMean = missingMean,
                     nonMissingMean = nonMissingMean,
                     pDiff = pDiff)    

    return(out)
}) %>% rbindlist 

missingnessMeans[, diffMeans := missingMean - nonMissingMean]

fwrite(missingnessMeans, file = "/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/missingnessMeans.csv", col.names = T)

phenoOut <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/processedData/scansUnadjustedFinal.csv")

## obtain skewness
skewness <-  apply(phenoOut[, ..pixels], 2, skewness, na.rm = T)

skewnessDT <- data.table(pixel = names(skewness),
                         skewness = skewness) %>%
   .[, c("y", "x") := tstrsplit(pixel, "_") ] 

## plot using geom_tile

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/plots/skewnessPixels.png", width = 800, height = 800)
ggplot(skewnessDT, aes(x = as.integer(x), y = as.integer(y), fill = skewness)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  scale_y_reverse() +
  theme_minimal() +
  labs(title = "Skewness of Retinal Thickness Phenotypes")
dev.off()

## randomly select 100 pixels and plot histograms of values for all pixels on one page

set.seed(123)
randPix <- sample(pixels, 100)

plotPix <- phenoOut[, ..randPix] %>%
  melt(variable.name = "pixel") %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins=20) +
  facet_wrap(~pixel, scales = "free") +
  theme_minimal() +
  ggtitle("Histograms of Random 100 pixels")

 png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/phenotypeCleaning/output/plots/histogramsRandomPixels.png", width = 1800, height = 1800) 
  print(plotPix)
dev.off()


getAgostinoP <- function(x) {

 x <- sort(x[complete.cases(x)])
    n <- length(x) %>% as.numeric
    s3 <- (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
    y <- s3 * sqrt((n + 1) * (n + 3)/(6 * (n - 2)))
    b2 <- 3 * (n * n + 27 * n - 70) * (n + 1) * (n + 3)/((n -
        2) * (n + 5) * (n + 7) * (n + 9))
    w <- sqrt(-1 + sqrt(2 * (b2 - 1)))
    d <- 1/sqrt(log(w))
    a <- sqrt(2/(w * w - 1))
    z <- d * log(y/a + sqrt((y/a)^2 + 1))
    pval <- 2*pnorm(z, lower.tail = FALSE)
    return(pval)
   } 