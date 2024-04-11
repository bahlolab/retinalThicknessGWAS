## 4 results files:
# output_h2 = "/vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/univariate_h2.txt"
# output_intercept = "/vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/univariate_ldscIntercept.txt"
# output_h2_stratified "/vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/stratified_h2.txt"
# output_intercept_stratified "/vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/stratified_ldscIntercept.txt"


## h2 files are {pixel}}\t{h2}\t{se_h2}
## intercept files are {pixel}}\t{intercept}\t{se_intercept}

## read in each
## split phenotype into x and y 
## plot h2 and intercept using geom tile


library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)

## univariate first....
h2 <- fread("/vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/univariate_h2.txt") %>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert = T)]

intercept <- fread("/vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/univariate_ldscIntercept.txt") %>%
    .[, c("y", "x") := tstrsplit(pixel, "_", type.convert = T)]


png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/ldsc/output/plots/univariate_h2.png", width = 600, height = 600)
ggplot(h2) +
  geom_tile(aes(y = y, x = x, fill = h2)) +
   scale_y_reverse() +
 theme_bw() +
  theme(legend.position = "bottom") + 
  labs(fill="heritability") 
dev.off()

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/ldsc/output/plots/univariate_intercept.png", width = 600, height = 600)
ggplot(intercept) +
  geom_tile(aes(y = y, x = x, fill = intercept)) +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "bottom") + 
  labs(fill="LDSC intercept")
dev.off()

## stratified next....
h2_stratified <- fread("/vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/stratified_h2.txt") %>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert = T)]

intercept_stratified <- fread("/vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/stratified_ldscIntercept.txt") %>%
  .[, c("y", "x") := tstrsplit(pixel, "_", type.convert = T)]

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/ldsc/output/plots/stratified_h2.png", width = 600, height = 600)
ggplot(h2_stratified) +
  geom_tile(aes(y = y, x = x, fill = h2)) +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "bottom") + 
  labs(fill="heritability")
dev.off()

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/ldsc/output/plots/stratified_intercept.png", width = 600, height = 600)
ggplot(intercept_stratified) +
  geom_tile(aes(y = y, x = x, fill = intercept)) +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "bottom") + 
  labs(fill="LDSC intercept")
dev.off()

