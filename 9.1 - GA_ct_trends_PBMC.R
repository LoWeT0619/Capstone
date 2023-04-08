# Library:
library(GEOquery)
library(dplyr)
library(ggplot2)
library(lumi)
library(EpiDISH)

# # Load PBMC data from NCBI:
# gset <- getGEO("GSE110828", GSEMatrix =TRUE, getGPL=FALSE)
# if (length(gset) > 1) idx <- grep("GPL13534", attr(gset, "names")) else idx <- 1
# gset <- gset[[idx]]
# 
# pd_450K <- pData(gset)
# beta_450K <- exprs(gset)

# Data:
load("/home/liuwent/06-New_Model_AddOn/00-Dataset/pd_450K_sub.RData")
load("/home/liuwent/06-New_Model_AddOn/00-Dataset/beta_450K.RData")

# Extract the PE and GA:
pd_PBMC <- pd_450K_sub[c('preeclampsia','gestational.age')]
colnames(pd_PBMC) <- c('Sample_Group','GA')
for(i in 1:nrow(pd_PBMC)){
  if(pd_PBMC$Sample_Group[i]==' Yes'){
    pd_PBMC$Sample_Group[i]='Disease'
  }else{pd_PBMC$Sample_Group[i]='Controls'}}
rownames(pd_PBMC) <- pd_450K_sub$geo_accession
# save(pd_PBMC, file='pd_PBMC.RData')
# load("/home/liuwent/04b-cell_type_deconvolution/pd_PBMC.RData")

# Manipulate cell type proportions and cell type deconvolution references:
## Extract the original PBMC data cell type proportions:
ctp_PBMC <- pd_450K_sub[,c(5:11)]
colnames(ctp_PBMC) <- c('Bcell','CD4T','CD8T','Gran','Mono','NK','nRBC')
rownames(ctp_PBMC) <- rownames(pd_PBMC)

## Omit Gran and re-calculate the cell type proportion in original PBMC data:
ctp_PBMC_nogran <- data.frame(matrix(nrow=nrow(ctp_PBMC), ncol=ncol(ctp_PBMC)-1))
colnames(ctp_PBMC_nogran) <- c('Bcell','CD4T','CD8T','Mono','NK','nRBC')
rownames(ctp_PBMC_nogran) <- rownames(ctp_PBMC)
for (i in 1:nrow(ctp_PBMC_nogran)){
  ctp_PBMC_nogran$Bcell[i] = round(ctp_PBMC$Bcell[i]/(1-ctp_PBMC$Gran[i]),6)
  ctp_PBMC_nogran$CD4T[i] = round(ctp_PBMC$CD4T[i]/(1-ctp_PBMC$Gran[i]),6)
  ctp_PBMC_nogran$CD8T[i] = round(ctp_PBMC$CD8T[i]/(1-ctp_PBMC$Gran[i]),6)
  ctp_PBMC_nogran$Mono[i] = round(ctp_PBMC$Mono[i]/(1-ctp_PBMC$Gran[i]),6)
  ctp_PBMC_nogran$NK[i] = round(ctp_PBMC$NK[i]/(1-ctp_PBMC$Gran[i]),6)
  ctp_PBMC_nogran$nRBC[i] = round(ctp_PBMC$nRBC[i]/(1-ctp_PBMC$Gran[i]),6)
}
# save(ctp_PBMC_nogran, file='ctp_PBMC_nogran.RData')
# load("/home/liuwent/04b-cell_type_deconvolution/ctp_PBMC_nogran.RData")

# ## PBMC cell type proportion by using newnewSetofMarkers (151,794) version:
# load("/home/liuwent/04b-cell_type_deconvolution/newnewSetofMarkers_final.RData")
# out_PBMC_150K<-epidish(beta.m=beta2m(beta_450K),ref.m=as.matrix(newnewSetofMarkers_final),method="CP")
# estF_PBMC_150K <- round(out_PBMC_150K$estF,6)
# 
# ## Omit Gran and re-calculate the PBMC cell type proportion by using newnewSetofMarkers (151,794) version:
# load("/home/liuwent/04b-cell_type_deconvolution/newnewSetofMarkers_final.RData")
# newnewSetofMarkers_final_nogran <- newnewSetofMarkers_final[, -6]
# out_PBMC_150K_nogran<-epidish(beta.m=beta2m(beta_450K),ref.m=as.matrix(newnewSetofMarkers_final_nogran),method="CP")
# estF_PBMC_150K_nogran <- round(out_PBMC_150K_nogran$estF,6)

# Scatterplot for GA and each cell types:
################################################################################
## Scatterplot for GA and each cell types (Grouped):
data <- data.frame(cbind(pd_PBMC$Sample_Group, pd_PBMC$GA, estF_PBMC_150K_nogran))
names(data)[names(data) == 'V1'] <- 'Sample_Group'
names(data)[names(data) == 'V2'] <- 'GA'

pdf("/home/liuwent/04b-cell_type_deconvolution/PBMC_GA_Bcell_scatterplot.pdf")
ggplot(data = data, aes(x = as.numeric(GA), y = as.numeric(Bcell), color=Sample_Group)) + 
  geom_point() +
  geom_smooth(aes(group=Sample_Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  labs(x = "GA", y = "Bcell")
dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/PBMC_GA_CD4T_scatterplot.pdf")
ggplot(data = data, aes(x = as.numeric(GA), y = as.numeric(CD4T), color=Sample_Group)) + 
  geom_point() +
  geom_smooth(aes(group=Sample_Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  labs(x = "GA", y = "CD4T")
dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/PBMC_GA_CD8T_scatterplot.pdf")
ggplot(data = data, aes(x = as.numeric(GA), y = as.numeric(CD8T), color=Sample_Group)) + 
  geom_point() +
  geom_smooth(aes(group=Sample_Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  labs(x = "GA", y = "CD8T")
dev.off()

# pdf("/home/liuwent/04b-cell_type_deconvolution/PBMC_GA_Gran_scatterplot.pdf")
# ggplot(data = data, aes(x = as.numeric(GA), y = as.numeric(Gran), color=Sample_Group)) +
#   geom_point() +
#   geom_smooth(aes(group=Sample_Group), method='lm') +
#   scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
#   labs(x = "GA", y = "Gran")
# dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/PBMC_GA_Mono_scatterplot.pdf")
ggplot(data = data, aes(x = as.numeric(GA), y = as.numeric(Mono), color=Sample_Group)) + 
  geom_point() +
  geom_smooth(aes(group=Sample_Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  labs(x = "GA", y = "Mono")
dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/PBMC_GA_NK_scatterplot.pdf")
ggplot(data = data, aes(x = as.numeric(GA), y = as.numeric(NK), color=Sample_Group)) + 
  geom_point() +
  geom_smooth(aes(group=Sample_Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  labs(x = "GA", y = "NK")
dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/PBMC_GA_nRBC_scatterplot.pdf")
ggplot(data = data, aes(x = as.numeric(GA), y = as.numeric(nRBC), color=Sample_Group)) + 
  geom_point() +
  geom_smooth(aes(group=Sample_Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  labs(x = "GA", y = "nRBC")
dev.off()

# ## Scatterplot for GA and each cell types (Ungrouped):
# data <- data.frame(cbind(pd$GA, estF))
# names(data)[names(data) == 'V1'] <- 'GA'
# 
# pdf("/home/liuwent/04b-cell_type_deconvolution/PBMC_GA_Bcell_scatterplot.pdf")
# ggplot(data = data, aes(x = GA, y = as.numeric(Bcell))) + 
#   geom_point() +
#   geom_smooth(method='lm') +
#   scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
#   labs(x = "GA", y = "Bcell")
# dev.off()
# 
# pdf("/home/liuwent/04b-cell_type_deconvolution/PBMC_GA_CD4T_scatterplot.pdf")
# ggplot(data = data, aes(x = GA, y = as.numeric(CD4T))) + 
#   geom_point() +
#   geom_smooth(method='lm') +
#   scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
#   labs(x = "GA", y = "CD4T")
# dev.off()
# 
# pdf("/home/liuwent/04b-cell_type_deconvolution/PBMC_GA_CD8T_scatterplot.pdf")
# ggplot(data = data, aes(x = GA, y = as.numeric(CD8T))) + 
#   geom_point() +
#   geom_smooth(method='lm') +
#   scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
#   labs(x = "GA", y = "CD8T")
# dev.off()
# 
# pdf("/home/liuwent/04b-cell_type_deconvolution/PBMC_GA_Gran_scatterplot.pdf")
# ggplot(data = data, aes(x = GA, y = as.numeric(Gran))) + 
#   geom_point() +
#   geom_smooth(method='lm') +
#   scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
#   labs(x = "GA", y = "Gran")
# dev.off()
# 
# pdf("/home/liuwent/04b-cell_type_deconvolution/PBMC_GA_Mono_scatterplot.pdf")
# ggplot(data = data, aes(x = GA, y = as.numeric(Mono))) + 
#   geom_point() +
#   geom_smooth(method='lm') +
#   scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
#   labs(x = "GA", y = "Mono")
# dev.off()
# 
# pdf("/home/liuwent/04b-cell_type_deconvolution/PBMC_GA_NK_scatterplot.pdf")
# ggplot(data = data, aes(x = GA, y = as.numeric(NK))) + 
#   geom_point() +
#   geom_smooth(method='lm') +
#   scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
#   labs(x = "GA", y = "NK")
# dev.off()
# 
# pdf("/home/liuwent/04b-cell_type_deconvolution/PBMC_GA_nRBC_scatterplot.pdf")
# ggplot(data = data, aes(x = GA, y = as.numeric(nRBC))) + 
#   geom_point() +
#   geom_smooth(method='lm') +
#   scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
#   labs(x = "GA", y = "nRBC")
# dev.off()