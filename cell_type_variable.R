#heatmap

library(dplyr)
library(limma)
library(tidyr)
library(ggplot2)
library(ggpubr)


pd_all <- read.csv("~/GitHub/CB_DNAm_PE/pd_all.csv")
estF <- pd_all[, 16:22]
# fit for all significant confounders except PE by using out$estF:
pd = pd_all%>%dplyr::select("Sample_Group", "GA", "BMI", "Eth2", "Parity", "Age")
pdnames = names(pd)
formstr <- paste0(pdnames, collapse = ' + ')
formstr <- paste0('~', formstr)
formstr <- as.formula(formstr)
design = model.matrix(formstr, data=pd)
cell_type_var = lmFit(t(estF), design)
cell_type_var = eBayes(cell_type_var)
cell_type_var_pval = data.frame(cell_type_var$p.value)

pval = round(cell_type_var$p.value, 3)
pval
coef = round(cell_type_var$coefficients, 3)
coef


write.csv(pval, file = "cell_type_var_pval.csv")
write.csv(coef, file = "cell_type_var_coef.csv")

library(ggpubr)
# PE residual
ResidualsMatrix <- residuals(cell_type_var, t(estF))


mean_mat = matrix(rowMeans(t(estF)-outer(coef[,2],ifelse(pd$Sample_Group == "Disease", 1,0))-ResidualsMatrix), 7, 62)

PE_adjust = ResidualsMatrix +  outer(coef[,2],ifelse(pd$Sample_Group == "Disease", 1,0))+
  mean_mat


PE_adjust_df = data.frame(pd$Sample_Group, t(PE_adjust))
PE_longer = pivot_longer(PE_adjust_df, 2:8, names_to = "cell_type")
ggplot(PE_longer, aes(x = cell_type, y = value, fill = pd.Sample_Group))+geom_boxplot()+
  stat_compare_means(label = "p.signif")

cell_type_var_pval$cell_type = rownames(cell_type_var_pval)
colnames(cell_type_var_pval)[2] = "Preeclampsia"
colnames(cell_type_var_pval)[5] = "Caucasian"
colnames(cell_type_var_pval)[6] = "Pacific Islander"

# Pvalue Heatmap ---------------------------------------------------
cell_pval_long = 
  cell_type_var_pval%>%select(-X.Intercept.)%>%
  pivot_longer(1:7, names_to = "variable", values_to = "pval")%>%
  mutate(pval = -log10(pval))

cell_pval_long$Stars = rep("", nrow(cell_pval_long))
for(i in 1:nrow(cell_pval_long)){
  if(cell_pval_long$pval[i]>-log10(0.001)){
    cell_pval_long$Stars[i] = "***"
  }else if(cell_pval_long$pval[i]>-log10(0.01)){
    cell_pval_long$Stars[i] = "**"
  }else if(cell_pval_long$pval[i]> -log10(0.05)){
    cell_pval_long$Stars[i] = "*"
  }
}

cell_pval_long$variable = 
  factor(cell_pval_long$variable, 
         levels = c("Preeclampsia", "GA", "BMI", "Caucasian", "Pacific Islander", "Parity", "Age"))
#relevel(cell_pval_long$variable)
ggplot(cell_pval_long, aes(x = variable, y = cell_type, fill = pval))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "#FF3300")+
  geom_text(aes(label = Stars), col = "#000000")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(fill=guide_legend(title="-log10(pval)"))
  


# PE and GA
pd2 = cbind(pd, estF)
pe_ga_long = pivot_longer(pd2, 7:13, names_to = "cell_type", values_to = "prop")
ggplot(pe_ga_long%>%filter(cell_type =="Gran"), aes(x = as.factor(GA), y=prop, col = Sample_Group))+
  geom_point()+
  geom_smooth(method = "lm", formula = "prop~GA")

# PE and BMI
ggplot(pe_ga_long%>%filter(cell_type == "nRBC"), aes(x = BMI, y = prop, col = BMI))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_gradient(low = "blue", high = "#FF3300")+
  annotate("text", label = "p=0.0067", x = 40, y = 0.12, col = "blue")

ggplot(pe_ga_long%>%filter(cell_type == "CD8T"), aes(x = BMI, y = prop, col = BMI))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_gradient(low = "blue", high = "#FF3300")+
  annotate("text", label = "p=0.002", x = 40, y = 0.08, col = "blue")

ggplot(pe_ga_long%>%filter(cell_type == "CD8T"), aes(x = as.factor("Pacific Islander"), y = prop, fill = BMI))+
  geom_boxplot()+
  #scale_color_gradient(low = "blue", high = "#FF3300")+
  annotate("text", label = "p=0.002", x = 0.5, y = 0.08, col = "blue")

# cell type vs PE barplot----------------------------------

pd_noPE = pd_all%>%dplyr::select("GA", "BMI", "Eth2", "Parity", "Age")
rownames(pd_noPE) <- pd$Sample_Name
#rownames(pd_noPE) <- pd$Sample_Name
pdnames = names(pd_noPE)
formstr <- paste0(pdnames, collapse = ' + ')
formstr <- paste0('~', formstr)
formstr <- as.formula(formstr)
design = model.matrix(formstr, data=pd_noPE)
fit_noPE = lmFit(t(estF), design)
fit_noPE = eBayes(fit_noPE)
save(fit_noPE, file="fit_noPE.RData")

# residual + intercept 
Residuals <- residuals(fit_noPE, t(estF))
colnames(Residuals) = pd_all$Sample_Name
Residuals = data.frame(t(Residuals))
Residuals$Sample_Group = pd_all$Sample_Group

residual_long = pivot_longer(Residuals, 1:7, names_to = "cell_type", values_to = "residual")
ggplot(residual_long, aes(x = cell_type, y = residual, fill = Sample_Group))+
  geom_boxplot()+
  stat_compare_means(label = "p.signif", method = "t.test")

intercept = fit_noPE$coefficients[, 1]
df_intercept = data.frame(t(intercept))
df_intercept = rbind(df_intercept[rep(1,nrow(pd_all)), ])
cell_residual = t(Residuals)+df_intercept
rownames(cell_residual) = pd_all$Sample_Name
cell_residual$Sample_Group = pd_all$Sample_Group

# ggplot
cell_residual_long = 
  pivot_longer(cell_residual, 1:7, names_to = "cell_type", values_to = "residual_intercept")
ggplot(cell_residual_long, aes(x = cell_type, y = residual_intercept, fill = Sample_Group))+
  geom_boxplot()+
  stat_compare_means(label = "p.signif", method = "t.test")
