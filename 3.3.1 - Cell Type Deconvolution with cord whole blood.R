# Library:
library(ChAMP)
library(shiny)
library(ggplot2)
library(parallel)
library(dplyr)
library(lumi)
library(minfi)
library(RColorBrewer)
library(MASS)
library(limma)
library(ggcorrplot)
library(carData)
library(car)
# library(AnalyzeFMRI)
library(EpiDISH)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annEPIC = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Shortcut
load("/home/liuwent/04-Full_Model/myLoad.Rdata")
load("/home/liuwent/04-Full_Model/pd.RData")
load("/home/liuwent/04-Full_Model/myCombat.RData")
load("/home/liuwent/04-Full_Model/Mvalues.RData")
load("/home/liuwent/04-Full_Model/CBref_lin.RData")

# Cell Type Deconvolution with cord whole blood (67 samples)
out<-epidish(beta.m=myCombat, ref.m=as.matrix(CBref),method="CP")
count <- out$estF
a = count[myLoad$pd$Sample_Group=="Controls",]
b = count[myLoad$pd$Sample_Group=="Disease",]
pdf("/home/liuwent/04-Full_Model/23-decon_linNRBC_CP.pdf")
boxplot(a, at=0:6*3 + 1, xlim=c(0, 21), ylim=range(a, b), xaxt="n", 
        col="red", main="", ylab="Cell type proportion")
boxplot(b, at=0:6*3 + 2, xaxt="n", add=TRUE, col="blue")
axis(1, at=0:6*3 + 1.5, labels=colnames(a), tick=TRUE)
text(x=7.5, y=0.1, "**")
text(x=10.5, y=0.12, "**")
text(x=16.5, y=max(count[,"Gran"]), "**")
text(x=19.5, y=0.2, "**")
legend("topleft", legend=c("Control_CP", "Cases_CP"), fill=c("red","blue"))
dev.off()

# Markers version:
library(EpiDISH)
library(ggplot2)
library(reshape2)

load("/home/liuwent/04-Full_Model/pd.RData")
load("/home/liuwent/04-Full_Model/myCombat.RData")
load("/home/liuwent/04b-cell_type_deconvolution/ct_ref_markers.RData")

out<-epidish(beta.m=myCombat, ref.m=as.matrix(ct_ref_markers),method="CP")
count <- out$estF

df = data.frame(cbind(count, pd))

Sample_Group1 = df$Sample_Group
df1 = data.frame(cbind(df[, 1:7], Sample_Group1))

df1_long = melt(df1, id.vars = "Sample_Group1")

pdf("/home/liuwent/04b-cell_type_deconvolution/23-ct_ref_markers_CP.pdf")
ggplot(df1_long, aes(variable, value, fill=Sample_Group1)) + 
  geom_boxplot() + 
  labs(title="Cell Type Proportion by Sample Group")
dev.off()

# newSetofMarkers version:
library(EpiDISH)
library(ggplot2)
library(reshape2)

load("/home/liuwent/04-Full_Model/pd.RData")
load("/home/liuwent/04-Full_Model/myCombat.RData")
load("/home/liuwent/04b-cell_type_deconvolution/newSetofMarkers_final.RData")

out<-epidish(beta.m=myCombat, ref.m=as.matrix(newSetofMarkers_final),method="CP")
count <- out$estF

df = data.frame(cbind(count, pd))

Sample_Group1 = df$Sample_Group
df1 = data.frame(cbind(df[, 1:7], Sample_Group1))

df1_long = melt(df1, id.vars = "Sample_Group1")

pdf("/home/liuwent/04b-cell_type_deconvolution/23-newSetofMarkers_CP.pdf")
ggplot(df1_long, aes(variable, value, fill=Sample_Group1)) + 
  geom_boxplot() + 
  labs(title="Cell Type Proportion by Sample Group")
dev.off()

# newSetofMarkers (868) version:
library(EpiDISH)
library(ggplot2)
library(reshape2)

load("/home/liuwent/04-Full_Model/pd.RData")
load("/home/liuwent/04-Full_Model/myCombat.RData")
load("/home/liuwent/04b-cell_type_deconvolution/newSetofMarkers_final.RData")

out<-epidish(beta.m=myCombat, ref.m=as.matrix(newSetofMarkers_final),method="CP")
count <- out$estF

df = data.frame(cbind(count, pd))

Sample_Group1 = df$Sample_Group
df1 = data.frame(cbind(df[, 1:7], Sample_Group1))

df1_long = melt(df1, id.vars = "Sample_Group1")

pdf("/home/liuwent/04b-cell_type_deconvolution/23-newSetofMarkers_CP(868).pdf")
ggplot(df1_long, aes(variable, value, fill=Sample_Group1)) + 
  geom_boxplot() + 
  labs(title="Cell Type Proportion by Sample Group")
dev.off()

# newSetofMarkers version:
library(EpiDISH)
library(ggplot2)
library(reshape2)

load("/home/liuwent/04-Full_Model/pd.RData")
load("/home/liuwent/04-Full_Model/myCombat.RData")
load("/home/liuwent/04b-cell_type_deconvolution/newSetofMarkers_final.RData")

out<-epidish(beta.m=myCombat, ref.m=as.matrix(newSetofMarkers_final),method="CP")
count <- out$estF

df = data.frame(cbind(count, pd))

Sample_Group1 = df$Sample_Group
df1 = data.frame(cbind(df[, 1:7], Sample_Group1))

df1_long = melt(df1, id.vars = "Sample_Group1")

pdf("/home/liuwent/04b-cell_type_deconvolution/23-newSetofMarkers_CP.pdf")
ggplot(df1_long, aes(variable, value, fill=Sample_Group1)) + 
  geom_boxplot() + 
  labs(title="Cell Type Proportion by Sample Group")
dev.off()

# newnewSetofMarkers (151,794) version:
library(EpiDISH)
library(ggplot2)
library(reshape2)
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
count <- out$estF

df = data.frame(cbind(count, pd))

Sample_Group1 = df$Sample_Group
df1 = data.frame(cbind(df[, 1:7], Sample_Group1))

df1_long = melt(df1, id.vars = "Sample_Group1")

pdf("/home/liuwent/04b-cell_type_deconvolution/23-newnewSetofMarkers_CP(151,794).pdf")
ggplot(df1_long, aes(variable, value, fill=Sample_Group1)) + 
  geom_boxplot() +
  labs(title="Cell Type Proportion by Sample Group") + 
  stat_compare_means(aes(label = after_stat(p.signif)))
dev.off()