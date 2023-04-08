# newnewSetofMarkers (151,794) residual version:
# Library:
library(EpiDISH)
library(ggplot2)
library(reshape2)
library(dplyr)
library(limma)
library(ggpubr)
library(ggsignif)

load("/home/liuwent/04-Full_Model/pd.RData")
load("/home/liuwent/04-Full_Model/myCombat.RData")
load("/home/liuwent/04b-cell_type_deconvolution/newnewSetofMarkers_final.RData")

# load("/home/liuwent/04b-cell_type_deconvolution/newnewSetofMarkers_pd.RData")
# load("/home/liuwent/04b-cell_type_deconvolution/newnewSetofMarkers_beta.RData")
# 
# Bcell = rowMeans(newnewSetofMarkers_beta[, which(newnewSetofMarkers_pd$CellType == "Bcell")])
# CD4T = rowMeans(newnewSetofMarkers_beta[, which(newnewSetofMarkers_pd$CellType == "CD4T")])
# CD8T = rowMeans(newnewSetofMarkers_beta[, which(newnewSetofMarkers_pd$CellType == "CD8T")])
# Gran = rowMeans(newnewSetofMarkers_beta[, which(newnewSetofMarkers_pd$CellType == "Gran")])
# Mono = rowMeans(newnewSetofMarkers_beta[, which(newnewSetofMarkers_pd$CellType == "Mono")])
# NK = rowMeans(newnewSetofMarkers_beta[, which(newnewSetofMarkers_pd$CellType == "NK")])
# nRBC = rowMeans(newnewSetofMarkers_beta[, which(newnewSetofMarkers_pd$CellType == "nRBC")])
# 
# # newnewSetofMarkers_final = cbind(Bcell, CD4T, CD8T, Gran, Mono, NK, nRBC)
# newnewSetofMarkers_final = cbind(CD8T, CD4T, NK, Bcell, Mono, Gran, nRBC)
# dim(newnewSetofMarkers_final)
# # 381,923 or 151,794   7
# # save(newnewSetofMarkers_final, file = "newnewSetofMarkers_final.RData")
# # load("/home/liuwent/04b-cell_type_deconvolution/newnewSetofMarkers_final.RData")

out<-epidish(beta.m=myCombat, ref.m=as.matrix(newnewSetofMarkers_final),method="CP")
estF <- out$estF

# fit for all significant confounders except PE by using out$estF:
pd_noPE = pd%>%dplyr::select("GA", "BMI", "Eth2", "Parity", "Age")
pdnames = names(pd_noPE)
formstr <- paste0(pdnames, collapse = ' + ')
formstr <- paste0('~', formstr)
formstr <- as.formula(formstr)
design = model.matrix(formstr, data=pd_noPE)
fit_noPE = lmFit(t(estF), design)
fit_noPE = eBayes(fit_noPE)
# save(fit_noPE, file="fit_noPE.RData")

Residuals_newnew_markers <- residuals(fit_noPE, t(estF))
# Add mean:
Residuals_newnew_markers <- Residuals_newnew_markers+matrix(apply(t(estF), 1, mean),
                                                            nrow=ncol(estF),
                                                            ncol=nrow(estF))
# # Add intercept:
# Residuals_newnew_markers <- Residuals_newnew_markers + fit_noPE$coefficients[,1]
# # save(Residuals_newnew_markers, file="Residuals_newnew_markers.RData")

count <- t(Residuals_newnew_markers)

df = data.frame(cbind(count, pd))

Sample_Group1 = df$Sample_Group
df1 = data.frame(cbind(df[, 1:7], Sample_Group1))

df1_long = melt(df1, id.vars = "Sample_Group1")

pdf("/home/liuwent/04b-cell_type_deconvolution/23-Residuals_newnewSetofMarkers_CP(151,794).pdf")
ggplot(df1_long, aes(variable, value, fill=Sample_Group1)) + 
  geom_boxplot() + 
  labs(title="Cell Type Residuals by Sample Group") + 
  stat_compare_means(aes(label = after_stat(p.signif)))
dev.off()