
library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)
library(patchwork)


origResults <- fread("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/finalResultsEUR/gwSigLociSummary.csv") %>%
.[BonferroniSig == "Y"]

# read in the results from the sensitivity analyses
# smoking results in /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/sensitivityAnalysesGWAS/output/chr${chr}Pixel${pixel}_smoking.${pixel}.glm.linear
# noSurgery results in /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/sensitivityAnalysesGWAS/output/chr${chr}Pixel${pixel}_noSurgery.${pixel}.glm.linear
# for each row in origResults, read in the corresponding sensitivity analysis results
# extract the sensitivity analysis results for that SNP 


smoking <- lapply(c(1:nrow(origResults)), function(locus) {

    print(paste0("Processing locus ", locus, " of ", nrow(origResults)))

    chr <- ifelse(origResults[locus, "CHR"] == 23, "X", origResults[locus, "CHR"])
    pixel <- origResults[locus, "pixel"]
    SNP <- origResults[locus, "ID"]

    file <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/sensitivityAnalysesGWAS/output/chr", chr, "Pixel", pixel, "_smoking.",pixel,".glm.linear")

    res <- fread(file) %>%
    .[ID==SNP] %>%
    .[, .(ID, A1, BETA, SE, P)] %>%
    setnames(., c("ID", "A1", "BETA_smoking", "SE_smoking", "P_smoking"))

    return(res)

}) %>%
rbindlist

noSurgery <- lapply(c(1:nrow(origResults)), function(locus) {

    print(paste0("Processing locus ", locus, " of ", nrow(origResults)))

    chr <- ifelse(origResults[locus, "CHR"] == 23, "X", origResults[locus, "CHR"])
    pixel <- origResults[locus, "pixel"]
    SNP <- origResults[locus, "ID"]

    file <- paste0("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/sensitivityAnalysesGWAS/output/chr", chr, "Pixel", pixel, "_noSurgery.",pixel,".glm.linear")

    res <- fread(file) %>%
    .[ID==SNP] %>%
    .[, .(ID, A1, BETA, SE, P)] %>%
    setnames(., c("ID", "A1", "BETA_noSurgery", "SE_noSurgery", "P_noSurgery"))

    return(res)

}) %>%
rbindlist

## merge origResults with smoking and noSurgery, on "ID" and "A1".
## plot comparison of betas and -log10 p-values for origResults vs smoking and origResults vs noSurgery

results <- merge(origResults, smoking, by = c("ID", "A1")) %>%
merge(., noSurgery, by = c("ID", "A1")) 

origSmokeBetas <- ggplot(results, aes(x = BETA, y = BETA_smoking)) +
  geom_point(size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(title = "Original vs Smoking Sensitivity Analysis Betas",
       x = "Original Beta",
       y = "Smoking Sensitivity Analysis Beta")

origSmokeP <- ggplot(results, aes(x = -log10(P), y = -log10(P_smoking))) +
  geom_point(size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(title = "Original vs Smoking Sensitivity Analysis -log10(p)",
       x = "Original -log10(p)",
       y = "Smoking Sensitivity Analysis -log10(p)")    

## plot together using patchwork
png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/sensitivityAnalysesGWAS/output/plots/origVsSmokeComparison.png", width = 1200, height = 600)
print(origSmokeBetas + origSmokeP)
dev.off()

## repeat for noSurgery

origNoSurgeryBetas <- ggplot(results, aes(x = BETA, y = BETA_noSurgery)) +
  geom_point(size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(title = "Original vs No Surgery Sensitivity Analysis Betas",
       x = "Original Beta",
       y = "No Surgery Sensitivity Analysis Beta")

origNoSurgeryP <- ggplot(results, aes(x = -log10(P), y = -log10(P_noSurgery))) +
  geom_point(size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    labs(title = "Original vs No Surgery Sensitivity Analysis -log10(p)",
         x = "Original -log10(p)",
         y = "No Surgery Sensitivity Analysis -log10(p)")

## plot together using patchwork
png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/sensitivityAnalysesGWAS/output/plots/origVsNoSurgeryComparison.png", width = 1200, height = 600)
print(origNoSurgeryBetas + origNoSurgeryP)
dev.off()
