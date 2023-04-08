# library:
library(EpiDISH)
library(ggplot2)
library(reshape2)
library(dplyr)
library(limma)

load("/home/liuwent/04-Full_Model/myLoad.Rdata")
load("/home/liuwent/04-Full_Model/pd.RData")
load("/home/liuwent/04-Full_Model/myCombat.RData")
load("/home/liuwent/04-Full_Model/Mvalues.RData")
load("/home/liuwent/04b-cell_type_deconvolution/order_new.Rdata")

pd = order_new
allConfounders = pd%>%dplyr::select("Sample_Group", "GA", "Age", "Parity", "Eth2", 
                                    "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC")

ind <- sapply(allConfounders, is.numeric)
#scale numerical variables
f = function(x){scale(x, center = FALSE)}
allConfounders[ind] <- lapply(allConfounders[ind],f)
formula = paste0(names(allConfounders), collapse = ' + ')
formula = paste0("~", formula)
formula = formula(formula)
#design matrix
design2 = model.matrix(formula, data = allConfounders)
#???#default setting control = 1, need to set case = 1???#
# design2[,2] = 1-design2[,2]
colnames(design2)[2] = "Cases"
#fit linear model
# fitConfounders = lmFit(Mvalues[, which(!is.na(allConfounders$BMI))], design2)
fitConfounders = lmFit(Mvalues, design2)
fitConfounders = eBayes(fitConfounders)
# save(fitConfounders, file = "fitConfounders.RData")

#make a toptable on second column sample_group
allg.limma <- topTable(fitConfounders, coef="Cases", n=dim(fitConfounders)[1], adjust.method="fdr")
sigg.limma <- subset(allg.limma, adj.P.Val < 0.05)
save(sigg.limma, file='sigg.limma.RData')
nonsigg.limma <- subset(allg.limma, adj.P.Val >= 0.05)
dim(sigg.limma)[1]

# Make a volcano plot (label DMP with different colors)
pval2 = p.adjust(fitConfounders$p.value[,2], method = "fdr")
col_list2 = sapply(seq(1, dim(fitConfounders$coefficients)[1]), function(i) {
  if (pval2[i]<0.05) {return("red")}
  else {return("black")} })
cpg <- data.frame(logFC=fitConfounders$coefficients[,2],pval=pval2)
sig_hyper <- rownames(cpg[cpg$logFC>0 & cpg$pval<0.05,])
length(sig_hyper)
# save(sig_hyper, file='sig_hyper.RData')
sig_hypo <- rownames(cpg[cpg$logFC<0 & cpg$pval<0.05,])
# save(sig_hypo, file='sig_hypo.RData')
length(sig_hypo)

pdf("/home/liuwent/04b-cell_type_deconvolution/06-fitConfounders_cvct.pdf", width = 10, height = 10)
volcanoplot(fitConfounders,coef=2,col=col_list2,highlight=5,names=rownames(fitConfounders$coefficients))
# abline(h=-log10(0.05), lty=2, col="red")
abline(v=-1, lty=2, col="grey")
abline(v=1, lty=2, col="grey")
text(-1, 5, "x=-1")
text(1, 5, "x=1")
dev.off()

length(pval2[pval2<0.05])
length(sig_hyper)
length(sig_hypo)