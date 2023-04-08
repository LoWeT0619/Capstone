# Library
library(dplyr)
library(limma)
library(ggplot2)
library(tidyverse)
library(ggpubr)

# Shortcut
load("/home/liuwent/04-Full_Model/pd.RData")
load("/home/liuwent/04-Full_Model/myCombat.RData")
load("/home/liuwent/04b-cell_type_deconvolution/estF.RData")
load("/home/liuwent/04b-cell_type_deconvolution/Residuals.RData")

# Omit Gran and re-calculate the cell type proportion by using newnewSetofMarkers (151,794) version:
estF <- as.data.frame(estF)
estF_nogran <- data.frame(matrix(nrow=nrow(estF), ncol=ncol(estF)-1))
colnames(estF_nogran) <- c('Bcell','CD4T','CD8T','Mono','NK','nRBC')
rownames(estF_nogran) <- rownames(estF)
for (i in 1:nrow(estF_nogran)){
  estF_nogran$Bcell[i] = round(estF$Bcell[i]/(1-estF$Gran[i]),6)
  estF_nogran$CD4T[i] = round(estF$CD4T[i]/(1-estF$Gran[i]),6)
  estF_nogran$CD8T[i] = round(estF$CD8T[i]/(1-estF$Gran[i]),6)
  estF_nogran$Mono[i] = round(estF$Mono[i]/(1-estF$Gran[i]),6)
  estF_nogran$NK[i] = round(estF$NK[i]/(1-estF$Gran[i]),6)
  estF_nogran$nRBC[i] = round(estF$nRBC[i]/(1-estF$Gran[i]),6)
}
# save(estF_nogran, file="estF_nogran.RData")
# load("/home/liuwent/04b-cell_type_deconvolution/estF_nogran.RData")

# Scatterplot for GA and each cell types:
################################################################################
## Scatterplot for GA and each cell types (Grouped):
data <- data.frame(cbind(pd$Sample_Group, pd$GA, estF))
names(data)[names(data) == 'pd.Sample_Group'] <- 'Sample_Group'
names(data)[names(data) == 'pd.GA'] <- 'GA'
data <- as.data.frame(data)

pdf("/home/liuwent/04b-cell_type_deconvolution/GA_Bcell_scatterplot.pdf")
ggplot(data = data, aes(x = as.numeric(GA), y = as.numeric(Bcell), color=Sample_Group)) + 
  geom_point() +
  geom_smooth(aes(group=Sample_Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  labs(x = "GA", y = "Bcell")
dev.off()
model_Bcell <- lm(data = data, Bcell~as.factor(data$Sample_Group)+as.numeric(data$GA)+Sample_Group*GA)
summary(model_Bcell)

pdf("/home/liuwent/04b-cell_type_deconvolution/GA_CD4T_scatterplot.pdf")
ggplot(data = data, aes(x = as.numeric(GA), y = as.numeric(CD4T), color=Sample_Group)) + 
  geom_point() +
  geom_smooth(aes(group=Sample_Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  labs(x = "GA", y = "CD4T")
dev.off()
model_CD4T <- lm(data = data, CD4T~as.factor(data$Sample_Group)+as.numeric(data$GA)+Sample_Group*GA)
summary(model_CD4T)

pdf("/home/liuwent/04b-cell_type_deconvolution/GA_CD8T_scatterplot.pdf")
ggplot(data = data, aes(x = as.numeric(GA), y = as.numeric(CD8T), color=Sample_Group)) + 
  geom_point() +
  geom_smooth(aes(group=Sample_Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  labs(x = "GA", y = "CD8T")
dev.off()
model_CD8T <- lm(data = data, CD8T~as.factor(data$Sample_Group)+as.numeric(data$GA)+Sample_Group*GA)
summary(model_CD8T)

pdf("/home/liuwent/04b-cell_type_deconvolution/GA_Gran_scatterplot.pdf")
ggplot(data = data, aes(x = as.numeric(GA), y = as.numeric(Gran), color=Sample_Group)) +
  geom_point() +
  geom_smooth(aes(group=Sample_Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  labs(x = "GA", y = "Gran")
dev.off()
model_Gran <- lm(data = data, Gran~as.factor(data$Sample_Group)+as.numeric(data$GA)+Sample_Group*GA)
summary(model_Gran)

pdf("/home/liuwent/04b-cell_type_deconvolution/GA_Mono_scatterplot.pdf")
ggplot(data = data, aes(x = as.numeric(GA), y = as.numeric(Mono), color=Sample_Group)) + 
  geom_point() +
  geom_smooth(aes(group=Sample_Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  labs(x = "GA", y = "Mono")
dev.off()
model_Mono <- lm(data = data, Mono~as.factor(data$Sample_Group)+as.numeric(data$GA)+Sample_Group*GA)
summary(model_Mono)

pdf("/home/liuwent/04b-cell_type_deconvolution/GA_NK_scatterplot.pdf")
ggplot(data = data, aes(x = as.numeric(GA), y = as.numeric(NK), color=Sample_Group)) + 
  geom_point() +
  geom_smooth(aes(group=Sample_Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  labs(x = "GA", y = "NK")
dev.off()
model_NK <- lm(data = data, NK~as.factor(data$Sample_Group)+as.numeric(data$GA)+Sample_Group*GA)
summary(model_NK)

pdf("/home/liuwent/04b-cell_type_deconvolution/GA_nRBC_scatterplot.pdf")
ggplot(data = data, aes(x = as.numeric(GA), y = as.numeric(nRBC), color=Sample_Group)) + 
  geom_point() +
  geom_smooth(aes(group=Sample_Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  labs(x = "GA", y = "nRBC")
dev.off()
model_nRBC <- lm(data = data, nRBC~as.factor(data$Sample_Group)+as.numeric(data$GA)+Sample_Group*GA)
summary(model_nRBC)

## Scatterplot for GA and each cell types (Ungrouped):
data <- data.frame(cbind(pd$GA, estF))
names(data)[names(data) == 'V1'] <- 'GA'

pdf("/home/liuwent/04b-cell_type_deconvolution/GA_Bcell_scatterplot.pdf")
ggplot(data = data, aes(x = GA, y = as.numeric(Bcell))) + 
  geom_point() +
  geom_smooth(method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  labs(x = "GA", y = "Bcell")
dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/GA_CD4T_scatterplot.pdf")
ggplot(data = data, aes(x = GA, y = as.numeric(CD4T))) + 
  geom_point() +
  geom_smooth(method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  labs(x = "GA", y = "CD4T")
dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/GA_CD8T_scatterplot.pdf")
ggplot(data = data, aes(x = GA, y = as.numeric(CD8T))) + 
  geom_point() +
  geom_smooth(method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  labs(x = "GA", y = "CD8T")
dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/GA_Gran_scatterplot.pdf")
ggplot(data = data, aes(x = GA, y = as.numeric(Gran))) + 
  geom_point() +
  geom_smooth(method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  labs(x = "GA", y = "Gran")
dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/GA_Mono_scatterplot.pdf")
ggplot(data = data, aes(x = GA, y = as.numeric(Mono))) + 
  geom_point() +
  geom_smooth(method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  labs(x = "GA", y = "Mono")
dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/GA_NK_scatterplot.pdf")
ggplot(data = data, aes(x = GA, y = as.numeric(NK))) + 
  geom_point() +
  geom_smooth(method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  labs(x = "GA", y = "NK")
dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/GA_nRBC_scatterplot.pdf")
ggplot(data = data, aes(x = GA, y = as.numeric(nRBC))) + 
  geom_point() +
  geom_smooth(method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  labs(x = "GA", y = "nRBC")
dev.off()

# Scatterplot for GA and each cell types (Residual Version):
################################################################################
data <- data.frame(cbind(pd$Sample_Group, pd$GA, t(Residuals)))
names(data)[names(data) == 'V1'] <- 'Sample_Group'
names(data)[names(data) == 'V2'] <- 'GA'

pdf("/home/liuwent/04b-cell_type_deconvolution/GA_Bcell_scatterplot.pdf")
ggplot(data = data, aes(x = GA, y = as.numeric(Bcell), color=Sample_Group)) + 
  geom_point() +
  geom_smooth(aes(group=Sample_Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  labs(x = "GA", y = "Bcell")
dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/GA_CD4T_scatterplot.pdf")
ggplot(data = data, aes(x = GA, y = as.numeric(CD4T), color=Sample_Group)) + 
  geom_point() +
  geom_smooth(aes(group=Sample_Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  labs(x = "GA", y = "CD4T")
dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/GA_CD8T_scatterplot.pdf")
ggplot(data = data, aes(x = GA, y = as.numeric(CD8T), color=Sample_Group)) + 
  geom_point() +
  geom_smooth(aes(group=Sample_Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  labs(x = "GA", y = "CD8T")
dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/GA_Gran_scatterplot.pdf")
ggplot(data = data, aes(x = GA, y = as.numeric(Gran), color=Sample_Group)) + 
  geom_point() +
  geom_smooth(aes(group=Sample_Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  labs(x = "GA", y = "Gran")
dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/GA_Mono_scatterplot.pdf")
ggplot(data = data, aes(x = GA, y = as.numeric(Mono), color=Sample_Group)) + 
  geom_point() +
  geom_smooth(aes(group=Sample_Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  labs(x = "GA", y = "Mono")
dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/GA_NK_scatterplot.pdf")
ggplot(data = data, aes(x = GA, y = as.numeric(NK), color=Sample_Group)) + 
  geom_point() +
  geom_smooth(aes(group=Sample_Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  labs(x = "GA", y = "NK")
dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/GA_nRBC_scatterplot.pdf")
ggplot(data = data, aes(x = GA, y = as.numeric(nRBC), color=Sample_Group)) + 
  geom_point() +
  geom_smooth(aes(group=Sample_Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  labs(x = "GA", y = "nRBC")
dev.off()