# Library:
library(ggplot2)
library(lumi)
library(limma)

# Shortcut
load("/home/liuwent/04-Full_Model/myLoad.Rdata")
load("/home/liuwent/04-Full_Model/pd.RData")
load("/home/liuwent/04-Full_Model/myCombat.RData")
load("/home/liuwent/04-Full_Model/Mvalues.RData")

## model2: regress only on sample group-----------------------------
myLoad$pd$Sample_Group = relevel(factor(myLoad$pd$Sample_Group), ref="Controls")
design0 = model.matrix(~myLoad$pd$Sample_Group)
fit0 = lmFit(beta2m(myCombat), design0)
fit00 = eBayes(fit0)
save(fit00, file = "fit00.RData")

## Get log_2 fold change and p value
log2fc = fit00$coefficients[,2]
##pval = p.adjust(fit$p.value[,2], method="bonferroni")
pval0 = p.adjust(fit00$p.value[,2], "BH")
length(pval0[pval0<0.05])
# 229,730

cpg <- data.frame(logFC=log2fc,pval=pval0)
sig_hyper <- rownames(cpg[cpg$logFC>0 & cpg$pval<0.05,])
length(sig_hyper)
# save(sig_hyper, file='sig_hyper.RData')
sig_hypo <- rownames(cpg[cpg$logFC<0 & cpg$pval<0.05,])
# save(sig_hypo, file='sig_hypo.RData')
length(sig_hypo)

table_fit00 = topTable(fit00, num=Inf, coef=2, adjust.method = "BH")
table_fit00$Type = ifelse(table_fit00$logFC>0, "Hyper", "Hypo")
length(table_fit00$adj.P.Val<0.05) # 24,353

## Make a volcano plot (label DMP with different colors)
col_list = sapply(seq(1, dim(fit00$coefficients)[1]), function(i) {
  if (pval0[i]<0.05) {return("red")}
  else {return("black")} })

pdf(file = "/home/liuwent/04-Full_Model/07-fit00.pdf", width = 10, height = 10)
volcanoplot(fit00, coef=2, col=col_list, highlight=5,
            names=rownames(fit00$coefficients))
# abline(h=-log10(0.05), lty=2, col="red")
abline(v=-1, lty=2, col="grey")
abline(v=1, lty=2, col="grey")
text(-1, 5, "x=-1")
text(1, 5, "x=1")
dev.off()