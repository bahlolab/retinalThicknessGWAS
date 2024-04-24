#!/usr/bin/env Rscript

library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)
library(here)

## univariate first....

smoking <- lapply(c(1:6), function(i) {
    fread(here("ldscOutFiles", "fPCs", paste0("fPC",i,".smoking.log")), skip = "p1", nrows=1) %>%
    .[, fPC := paste0("fPC",i)]

}) %>%
rbindlist %>%
  .[, fPC := factor(fPC, levels = paste0("fPC", 6:1))]

noSurgery <- lapply(c(1:6), function(i) {
    fread(here("ldscOutFiles", "fPCs", paste0("fPC",i,".noSurgery.log")), skip = "p1", nrows=1) %>%
    .[, fPC := paste0("fPC",i)]

}) %>%
rbindlist %>%
  .[, fPC := factor(fPC, levels = paste0("fPC", 6:1))]


png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/ldscSensitivity/output/plots/smokingGenCorr.png", width = 600, height = 600)
ggplot(smoking) +
  geom_point(aes(y = fPC, x = rg)) +
  geom_errorbar(aes(y = fPC, xmin = rg - se, xmax = rg + se), width = 0.1) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  theme_bw() +
  scale_y_discrete(labels = levels(smoking$fPC)) +
  scale_x_continuous(limits = c(0.95, 1.01)) +
  labs(y = "fPC", x = "Genetic Correlation", title = "Smoking Sensitivity Analysis" )
dev.off()

png("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/ldscSensitivity/output/plots/noSurgeryGenCorr.png", width = 600, height = 600)
ggplot(noSurgery) +
  geom_point(aes(y = fPC, x = rg)) +
  geom_errorbar(aes(y = fPC, xmin = rg - se, xmax = rg + se), width = 0.1) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  theme_bw() +
  scale_y_discrete(labels = levels(noSurgery$fPC)) +
  scale_x_continuous(limits = c(0.95, 1.01)) +
  labs(y = "fPC", x = "Genetic Correlation", title = "Surgery exclusions Sensitivity Analysis" )
dev.off()

fwrite(smoking, here("ldscOutFiles", "fPCs", "smokingGenCorr.txt"))
fwrite(noSurgery, here("ldscOutFiles", "fPCs", "noSurgeryGenCorr.txt"))