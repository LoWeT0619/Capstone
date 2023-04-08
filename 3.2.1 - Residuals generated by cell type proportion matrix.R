# Library
library(ggplot2)
library(dplyr)
library(limma)
library(EpiDISH)

# Shortcut
load("/home/liuwent/04-Full_Model/myLoad.Rdata")
load("/home/liuwent/04-Full_Model/pd.RData")
load("/home/liuwent/04-Full_Model/myCombat.RData")
load("/home/liuwent/04-Full_Model/Mvalues.RData")
load("/home/liuwent/04b-cell_type_deconvolution/newnewSetofMarkers_final.RData")

out<-epidish(beta.m=myCombat, ref.m=as.matrix(newnewSetofMarkers_final),method="CP")
estF <- out$estF
save(estF, file="estF.RData")

# fit for all significant confounders except PE by using out$estF:
pd_noPE = pd%>%dplyr::select("GA", "BMI", "Eth2", "Parity", "Age")
rownames(pd) <- pd$Sample_Name
rownames(pd_noPE) <- pd$Sample_Name
pdnames = names(pd_noPE)
formstr <- paste0(pdnames, collapse = ' + ')
formstr <- paste0('~', formstr)
formstr <- as.formula(formstr)
design = model.matrix(formstr, data=pd_noPE)
fit_noPE = lmFit(t(estF), design)
fit_noPE = eBayes(fit_noPE)
save(fit_noPE, file="fit_noPE.RData")

Residuals <- residuals(fit_noPE, t(estF))
Residuals <- Residuals+matrix(apply(t(estF), 1, mean), nrow=ncol(estF), ncol=nrow(estF))
save(Residuals, file="Residuals.RData")
load("/home/liuwent/04b-cell_type_deconvolution/Residuals.RData")