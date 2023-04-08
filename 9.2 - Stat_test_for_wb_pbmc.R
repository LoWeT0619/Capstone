# Library:
library(lumi)
library(EpiDISH)
library(dplyr)
library(ggplot2)
library(limma)
library(tidyverse)

# Data:
load("/home/liuwent/04-Full_Model/pd.RData")
load("/home/liuwent/04b-cell_type_deconvolution/pd_PBMC.RData")

load("/home/liuwent/04b-cell_type_deconvolution/estF_nogran.RData") # whole blood (without gran)
data_wb <- data.frame(cbind(pd$Sample_Group, pd$GA, estF_nogran))
names(data_wb)[names(data_wb) == 'pd.Sample_Group'] <- 'Sample_Group'
names(data_wb)[names(data_wb) == 'pd.GA'] <- 'GA'
data_wb$Group <- 'WB'
data_wb_PE <- data_wb[which(data_wb$Sample_Group=='Disease'), ]
data_wb_Controls <- data_wb[which(data_wb$Sample_Group=='Controls'), ]

load("/home/liuwent/04b-cell_type_deconvolution/ctp_PBMC_nogran.RData") # PBMC (without gran)
data_pbmc <- data.frame(cbind(pd_PBMC$Sample_Group, pd_PBMC$GA, ctp_PBMC_nogran))
names(data_pbmc)[names(data_pbmc) == 'pd_PBMC.Sample_Group'] <- 'Sample_Group'
names(data_pbmc)[names(data_pbmc) == 'pd_PBMC.GA'] <- 'GA'
data_pbmc$Group <- 'PBMC'
data_pbmc_PE <- data_pbmc[which(data_pbmc$Sample_Group=='Disease'), ]
data_pbmc_Controls <- data_pbmc[which(data_pbmc$Sample_Group=='Controls'), ]

data <- rbind(data_wb, data_pbmc)
data_PE <- rbind(data_wb_PE, data_pbmc_PE)
data_Controls <- rbind(data_wb_Controls, data_pbmc_Controls)

pdf("/home/liuwent/04b-cell_type_deconvolution/WB_PBMC_GA_Bcell_scatterplot.pdf")
ggplot(data = data_Controls, aes(x = as.numeric(GA), y = as.numeric(Bcell), color=Group)) + 
  geom_point() +
  geom_smooth(aes(group=Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) + 
  scale_color_manual(values = c("#d95f02", "#7570b3"))
  labs(x = "GA", y = "Bcell")
dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/WB_PBMC_GA_CD4T_scatterplot.pdf")
ggplot(data = data_Controls, aes(x = as.numeric(GA), y = as.numeric(CD4T), color=Group)) + 
  geom_point() +
  geom_smooth(aes(group=Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  scale_color_manual(values = c("#d95f02", "#7570b3"))
  labs(x = "GA", y = "CD4T")
dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/WB_PBMC_GA_CD8T_scatterplot.pdf")
ggplot(data = data_Controls, aes(x = as.numeric(GA), y = as.numeric(CD8T), color=Group)) + 
  geom_point() +
  geom_smooth(aes(group=Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  scale_color_manual(values = c("#d95f02", "#7570b3"))
  labs(x = "GA", y = "CD8T")
dev.off()

# pdf("/home/liuwent/04b-cell_type_deconvolution/WB_PBMC_GA_Gran_scatterplot.pdf")
# ggplot(data = data_Controls, aes(x = as.numeric(GA), y = as.numeric(Gran), color=Group)) +
#   geom_point() +
#   geom_smooth(aes(group=Group), method='lm') +
#   scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
#   scale_color_manual(values = c("#d95f02", "#7570b3"))
#   labs(x = "GA", y = "Gran")
# dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/WB_PBMC_GA_Mono_scatterplot.pdf")
ggplot(data = data_Controls, aes(x = as.numeric(GA), y = as.numeric(Mono), color=Group)) + 
  geom_point() +
  geom_smooth(aes(group=Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  scale_color_manual(values = c("#d95f02", "#7570b3"))
  labs(x = "GA", y = "Mono")
dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/WB_PBMC_GA_NK_scatterplot.pdf")
ggplot(data = data_Controls, aes(x = as.numeric(GA), y = as.numeric(NK), color=Group)) + 
  geom_point() +
  geom_smooth(aes(group=Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  scale_color_manual(values = c("#d95f02", "#7570b3"))
  labs(x = "GA", y = "NK")
dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/WB_PBMC_GA_nRBC_scatterplot.pdf")
ggplot(data = data_Controls, aes(x = as.numeric(GA), y = as.numeric(nRBC), color=Group)) + 
  geom_point() +
  geom_smooth(aes(group=Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  scale_color_manual(values = c("#d95f02", "#7570b3"))
  labs(x = "GA", y = "nRBC")
dev.off()

# Statistical test for slopes:
## For Bcell:
### Original:
model_Bcell <- lm(data = data, Bcell~as.factor(data$Group)+as.numeric(data$GA))
summary(model_Bcell)
### Disease group only:
model_PE_Bcell <- lm(data = data_PE, Bcell~as.factor(data_PE$Group)+as.numeric(data_PE$GA)+GA*Group)
summary(model_PE_Bcell)
### Control group only:
model_PE_Bcell <- lm(data = data_Controls, Bcell~as.factor(data_Controls$Group)+as.numeric(data_Controls$GA)+GA*Group)
summary(model_PE_Bcell)

## For CD4T:
model_CD4T <- lm(data = data, CD4T~as.factor(data$Group)+as.numeric(data$GA))
summary(model_CD4T)
### Disease group only:
model_PE_CD4T <- lm(data = data_PE, CD4T~as.factor(data_PE$Group)+as.numeric(data_PE$GA)+GA*Group)
summary(model_PE_CD4T)
### Control group only:
model_PE_CD4T <- lm(data = data_Controls, CD4T~as.factor(data_Controls$Group)+as.numeric(data_Controls$GA)+GA*Group)
summary(model_PE_CD4T)

## For CD8T:
model_CD8T <- lm(data = data, CD8T~as.factor(data$Group)+as.numeric(data$GA))
summary(model_CD8T)
### Disease group only:
model_PE_CD8T <- lm(data = data_PE, CD8T~as.factor(data_PE$Group)+as.numeric(data_PE$GA)+GA*Group)
summary(model_PE_CD8T)
### Control group only:
model_PE_CD8T <- lm(data = data_Controls, CD8T~as.factor(data_Controls$Group)+as.numeric(data_Controls$GA)+GA*Group)
summary(model_PE_CD8T)

## For Mono:
model_Mono <- lm(data = data, Mono~as.factor(data$Group)+as.numeric(data$GA))
summary(model_Mono)
### Disease group only:
model_PE_Mono <- lm(data = data_PE, Mono~as.factor(data_PE$Group)+as.numeric(data_PE$GA)+GA*Group)
summary(model_PE_Mono)
### Control group only:
model_PE_Mono <- lm(data = data_Controls, Mono~as.factor(data_Controls$Group)+as.numeric(data_Controls$GA)+GA*Group)
summary(model_PE_Mono)

## For NK:
model_NK <- lm(data = data, NK~as.factor(data$Group)+as.numeric(data$GA))
summary(model_NK)
### Disease group only:
model_PE_NK <- lm(data = data_PE, NK~as.factor(data_PE$Group)+as.numeric(data_PE$GA)+GA*Group)
summary(model_PE_NK)
### Control group only:
model_PE_NK <- lm(data = data_Controls, NK~as.factor(data_Controls$Group)+as.numeric(data_Controls$GA)+GA*Group)
summary(model_PE_NK)

## For nRBC:
model_nRBC <- lm(data = data, nRBC~as.factor(data$Group)+as.numeric(data$GA)+GA*Group)
summary(model_nRBC)
### Disease group only:
model_PE_nRBC <- lm(data = data_PE, nRBC~as.factor(data_PE$Group)+as.numeric(data_PE$GA)+GA*Group)
summary(model_PE_nRBC)
### Control group only:
model_PE_nRBC <- lm(data = data_Controls, nRBC~as.factor(data_Controls$Group)+as.numeric(data_Controls$GA)+GA*Group)
summary(model_PE_nRBC)
