#!/usr/bin/env Rscript

library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)

## 4 results files:
# output_h2 = "/vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/fPCs_univariate_h2.txt"
# output_intercept = "/vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/fPCs_univariate_ldscIntercept.txt"
# output_h2_stratified "/vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/fPCs_stratified_h2.txt"
# output_intercept_stratified "/vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/fPCs_stratified_ldscIntercept.txt"


## h2 files are {fPC}}\t{h2}\t{se}
## intercept files are {fPC}}\t{intercept}\t{se}

## read in each
## plot h2 and intercept +/- se using geom point and geom errorbar
## make y labels fPC1 to fPC6 (top to bottom)


## univariate first....

h2 <- fread("/vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/fPCs_univariate_h2.txt") %>%
  .[, fPC := factor(fPC, levels = paste0("fPC", 6:1))]
intercept <- fread("/vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/fPCs_univariate_ldscIntercept.txt") %>%
  .[, fPC := factor(fPC, levels = paste0("fPC", 6:1))]

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/ldsc/output/plots/fPCsUnivariate_h2.png", width = 600, height = 600)
ggplot(h2) +
  geom_point(aes(y = fPC, x = h2)) +
  geom_errorbar(aes(y = fPC, xmin = h2 - se, xmax = h2 + se), width = 0.1) +
  theme_bw() +
  scale_y_discrete(labels = levels(h2$fPC)) +
  labs(y = "fPC", x = "heritability") 
dev.off()

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/ldsc/output/plots/fPCsUnivariate_intercept.png", width = 600, height = 600)
ggplot(intercept) +
  geom_point(aes(y = fPC, x = intercept)) +
  geom_errorbar(aes(y = fPC, xmin = intercept - se, xmax = intercept + se), width = 0.1) +
  theme_bw() +
  scale_y_discrete(labels = levels(h2$fPC)) +
  labs(y = "fPC", x = "intercept")
dev.off()

## stratified next....

h2_stratified <- fread("/vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/fPCs_stratified_h2.txt") %>%
  .[, fPC := factor(fPC, levels = paste0("fPC", 6:1))]
intercept_stratified <- fread("/vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/fPCs_stratified_ldscIntercept.txt") %>%
  .[, fPC := factor(fPC, levels = paste0("fPC", 6:1))]

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/ldsc/output/plots/fPCsStratified_h2.png", width = 600, height = 600)
ggplot(h2_stratified) +
  geom_point(aes(y = fPC, x = h2)) +
  geom_errorbar(aes(y = fPC, xmin = h2 - se, xmax = h2 + se), width = 0.1) +
  theme_bw() +
  scale_y_discrete(labels = levels(h2$fPC)) +
  labs(y = "fPC", x = "heritability")
dev.off()

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/ldsc/output/plots/fPCsStratified_intercept.png", width = 600, height = 600)
ggplot(intercept_stratified) +
  geom_point(aes(y = fPC, x = intercept)) +
  geom_errorbar(aes(y = fPC, xmin = intercept - se, xmax = intercept + se), width = 0.1) +
  theme_bw() +
  scale_y_discrete(labels = levels(h2$fPC)) +
  labs(y = "fPC", x = "intercept")
dev.off()

