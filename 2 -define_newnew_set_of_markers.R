# Library:
library(wateRmelon)
library(RPMM)
library(minfi)
library(RColorBrewer)

# Shortcut:
load("/home/liuwent/04b-cell_type_deconvolution/pd_Lin.RData")
load("/home/liuwent/04b-cell_type_deconvolution/pd_B.RData")

load("/home/liuwent/04b-cell_type_deconvolution/MethylSets.RData")
load("/home/liuwent/04b-cell_type_deconvolution/MethylSets_scale.RData")

load("/home/liuwent/04b-cell_type_deconvolution/comb_pd.RData")
load("/home/liuwent/04b-cell_type_deconvolution/comb_beta.RData")

load("/home/liuwent/04b-cell_type_deconvolution/newnewSetofMarkers_pd.RData")
load("/home/liuwent/04b-cell_type_deconvolution/newnewSetofMarkers_beta.RData")

# Dataset:
## two references' raw data:
### Lin:
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
load("/home/liuwent/01-data/FlowSorted.CordTissueAndBlood.EPIC.rda")
Lin = updateObject(FlowSorted.CordTissueAndBlood.EPIC)
Mset.Noob.Lin = preprocessNoob(Lin)   # class: MethylSet
msetMapToGenome_Lin = mapToGenome(Mset.Noob.Lin)   # class: GenomicMethylSet
msetGenoRatio_Lin = ratioConvert(msetMapToGenome_Lin)
pd_Lin = pData(msetGenoRatio_Lin)
beta_Lin = getBeta(msetMapToGenome_Lin)

### Bukulski:
library("FlowSorted.CordBlood.450k")
data(FlowSorted.CordBlood.450k)
B = FlowSorted.CordBlood.450k
Mset.Noob.B = preprocessNoob(B)   # class: MethylSet
msetMapToGenome_B = mapToGenome(Mset.Noob.B)   # class: GenomicMethylSet
msetGenoRatio_B  = ratioConvert(msetMapToGenome_B)
pd_B = pData(msetGenoRatio_B)
beta_B = getBeta(msetMapToGenome_B)

# step 1: use the common CpGs in Lin and B data for downstream.
## use combineArrays from minfi to combine two MethylSets:
MethylSets <- combineArrays(Mset.Noob.Lin, Mset.Noob.B)
dim(MethylSets)   #453093    226
# save(MethylSets, file="MethylSets.RData")
# load("/home/liuwent/04b-cell_type_deconvolution/MethylSets.RData")

## PCA check:
### Lin:
load("/home/liuwent/04b-cell_type_deconvolution/pd_Lin.RData")
pdSub = pd_Lin[which(pd_Lin$CellType!="Endothelial" &
                       pd_Lin$CellType!="Stromal" &
                       pd_Lin$CellType!="Epithelial"),]
MethylSets_Lin <- MethylSets[, which(colnames(MethylSets)%in%rownames(pdSub))]
msetsMapToGenome_Lin = mapToGenome(MethylSets_Lin)
MethylSets_Lin_beta = getBeta(msetsMapToGenome_Lin)

beta = t(MethylSets_Lin_beta)
beta_stand = sweep(beta,2,colMeans(beta),"-")
SVD <- svd(beta_stand)
svd_plot = cbind(as.data.frame(rownames(pdSub)),
                 as.data.frame(SVD$u[, 1]), 
                 as.data.frame(SVD$u[, 2]), 
                 as.data.frame(pdSub$CellType))
colnames(svd_plot) = c("ID","PC1","PC2","CellType")
col.pat = brewer.pal(n = length(unique(pdSub$CellType)), name = "Set2")

pdf("/home/liuwent/04b-cell_type_deconvolution/MethylSets_Lin_PCA.pdf")
ggplot(svd_plot, aes(x=PC1, y=PC2, color=CellType)) +
  scale_color_manual(values=col.pat)+
  geom_point()+
  theme_classic(base_size = 18)
dev.off()

### Bakulski:
load("/home/liuwent/04b-cell_type_deconvolution/pd_B.RData")
MethylSets_B <- MethylSets[, which(colnames(MethylSets)%in%rownames(pd_B))]
msetsMapToGenome_B = mapToGenome(MethylSets_B)
MethylSets_B_beta = getBeta(msetsMapToGenome_B)

beta = t(MethylSets_B_beta)
beta_stand = sweep(beta,2,colMeans(beta),"-")
SVD <- svd(beta_stand)
svd_plot = cbind(as.data.frame(rownames(pd_B)),
                 as.data.frame(SVD$u[, 1]), 
                 as.data.frame(SVD$u[, 2]), 
                 as.data.frame(pd_B$CellType))
colnames(svd_plot) = c("ID","PC1","PC2","CellType")
col.pat = brewer.pal(n = length(unique(pd_B$CellType)), name = "Set2")

pdf("/home/liuwent/04b-cell_type_deconvolution/MethylSets_B_PCA.pdf")
ggplot(svd_plot, aes(x=PC1, y=PC2, color=CellType)) +
  scale_color_manual(values=col.pat)+
  geom_point()+
  theme_classic(base_size = 18)
dev.off()

### Comb:
load("/home/liuwent/04b-cell_type_deconvolution/comb_pd.RData")
load("/home/liuwent/04b-cell_type_deconvolution/MethylSets.RData")
MethylSets <- MethylSets[, which(colnames(MethylSets)%in%rownames(comb_pd))]
msetsMapToGenome = mapToGenome(MethylSets)
MethylSets_beta = getBeta(msetsMapToGenome)

beta = t(MethylSets_beta)
beta_stand = sweep(beta,2,colMeans(beta),"-")
SVD <- svd(beta_stand)
svd_plot = cbind(as.data.frame(rownames(comb_pd)),
                 as.data.frame(SVD$u[, 1]), 
                 as.data.frame(SVD$u[, 2]), 
                 as.data.frame(comb_pd$CellType))
colnames(svd_plot) = c("ID","PC1","PC2","CellType")
col.pat = brewer.pal(n = length(unique(comb_pd$CellType)), name = "Set2")

pdf("/home/liuwent/04b-cell_type_deconvolution/MethylSets_comb_PCA.pdf")
ggplot(svd_plot, aes(x=PC1, y=PC2, color=CellType)) +
  scale_color_manual(values=col.pat)+
  geom_point()+
  theme_classic(base_size = 18)
dev.off()

#step 2: rescale both datasets so they have the same scales:
## use BMIQ from wateRmelon to rescale the MethylSets:
MethylSets_scale <- BMIQ(MethylSets)
dim(MethylSets_scale)   #453093    226
# save(MethylSets_scale, file="MethylSets_scale.RData")
# load("/home/liuwent/04b-cell_type_deconvolution/MethylSets_scale.RData")

## PCA check:
### Lin:
load("/home/liuwent/04b-cell_type_deconvolution/pd_Lin.RData")
pdSub = pd_Lin[which(pd_Lin$CellType!="Endothelial" &
                       pd_Lin$CellType!="Stromal" &
                       pd_Lin$CellType!="Epithelial"),]
MethylSets_scale_Lin <- MethylSets_scale[, which(colnames(MethylSets_scale)%in%rownames(pdSub))]

beta = t(MethylSets_scale_Lin)
beta_stand = sweep(beta,2,colMeans(beta),"-")
SVD <- svd(beta_stand)
svd_plot = cbind(as.data.frame(rownames(pdSub)),
                 as.data.frame(SVD$u[, 1]), 
                 as.data.frame(SVD$u[, 2]), 
                 as.data.frame(pdSub$CellType))
colnames(svd_plot) = c("ID","PC1","PC2","CellType")
col.pat = brewer.pal(n = length(unique(pdSub$CellType)), name = "Set2")

pdf("/home/liuwent/04b-cell_type_deconvolution/MethylSets_scale_Lin_PCA.pdf")
ggplot(svd_plot, aes(x=PC1, y=PC2, color=CellType)) +
  scale_color_manual(values=col.pat)+
  geom_point()+
  theme_classic(base_size = 18)
dev.off()

### Bakulski:
load("/home/liuwent/04b-cell_type_deconvolution/pd_B.RData")
MethylSets_scale_B <- MethylSets_scale[, which(colnames(MethylSets_scale)%in%rownames(pd_B))]

beta = t(MethylSets_scale_B)
beta_stand = sweep(beta,2,colMeans(beta),"-")
SVD <- svd(beta_stand)
svd_plot = cbind(as.data.frame(rownames(pd_B)),
                 as.data.frame(SVD$u[, 1]), 
                 as.data.frame(SVD$u[, 2]), 
                 as.data.frame(pd_B$CellType))
colnames(svd_plot) = c("ID","PC1","PC2","CellType")
col.pat = brewer.pal(n = length(unique(pd_B$CellType)), name = "Set2")

pdf("/home/liuwent/04b-cell_type_deconvolution/MethylSets_scale_B_PCA.pdf")
ggplot(svd_plot, aes(x=PC1, y=PC2, color=CellType)) +
  scale_color_manual(values=col.pat)+
  geom_point()+
  theme_classic(base_size = 18)
dev.off()

### Comb:
load("/home/liuwent/04b-cell_type_deconvolution/comb_pd.RData")
load("/home/liuwent/04b-cell_type_deconvolution/MethylSets_scale.RData")
MethylSets_scale <- MethylSets_scale[, which(colnames(MethylSets_scale)%in%rownames(comb_pd))]

beta = t(MethylSets_scale)
beta_stand = sweep(beta,2,colMeans(beta),"-")
SVD <- svd(beta_stand)
svd_plot = cbind(as.data.frame(rownames(comb_pd)),
                 as.data.frame(SVD$u[, 1]), 
                 as.data.frame(SVD$u[, 2]), 
                 as.data.frame(comb_pd$CellType))
colnames(svd_plot) = c("ID","PC1","PC2","CellType")
col.pat = brewer.pal(n = length(unique(comb_pd$CellType)), name = "Set2")

pdf("/home/liuwent/04b-cell_type_deconvolution/MethylSets_scale_omb_PCA.pdf")
ggplot(svd_plot, aes(x=PC1, y=PC2, color=CellType)) +
  scale_color_manual(values=col.pat)+
  geom_point()+
  theme_classic(base_size = 18)
dev.off()

#step 3: combine lin's cell types with nRBC from B-data
load("/home/liuwent/04b-cell_type_deconvolution/pd_Lin.RData")
pdSub = pd_Lin[which(pd_Lin$CellType!="Endothelial" &
                       pd_Lin$CellType!="Stromal" &
                       pd_Lin$CellType!="Epithelial"),]
pdsub = pdSub[,c("CellType","gender")]

load("/home/liuwent/04b-cell_type_deconvolution/pd_B.RData")
nRBC=pd_B[which(pd_B$CellType=="nRBC"),]
nRBC = nRBC[,c("CellType","Sex")]
colnames(nRBC)[2]="gender"

comb_pd = rbind(pdsub,nRBC)
# save(comb_pd, file="comb_pd.RData")
# load("/home/liuwent/04b-cell_type_deconvolution/comb_pd.RData")

comb_beta = MethylSets_scale[, which(colnames(MethylSets_scale)%in%rownames(comb_pd))]
# save(comb_beta, file="comb_beta.RData")
# load("/home/liuwent/04b-cell_type_deconvolution/comb_beta.RData")

## PCA check:
beta = t(comb_beta)
beta_stand = sweep(beta,2,colMeans(beta),"-")
SVD <- svd(beta_stand)
svd_plot = cbind(as.data.frame(rownames(comb_pd)),
                 as.data.frame(SVD$u[, 1]), 
                 as.data.frame(SVD$u[, 2]), 
                 as.data.frame(comb_pd$CellType))
colnames(svd_plot) = c("ID","PC1","PC2","CellType")
col.pat = brewer.pal(n = length(unique(comb_pd$CellType)), name = "Set2")

pdf("/home/liuwent/04b-cell_type_deconvolution/MethylSets_scale_comb_PCA.pdf")
ggplot(svd_plot, aes(x=PC1, y=PC2, color=CellType)) +
  scale_color_manual(values=col.pat)+
  geom_point()+
  theme_classic(base_size = 18)
dev.off()

# Define the function:
deconv_pair <- function(beta,pd,celltypes,p_val,n_prob=100){
  library(limma)
  library(lumi)
  splitit <- function(x) {
    split(seq_along(x), x)
  }
  newnewSetofMarkers <- list()
  celltype_index = splitit(pd$CellType)
  for(i in celltypes){
    sig_limma_list = list()
    for(j in celltypes){
      if(!i==j){
        pd$status <- NA
        pd[pd$CellType==i,]$status <- rep(0,dim(pd[pd$CellType==i,])[1])
        pd[pd$CellType==j,]$status <- rep(1,dim(pd[pd$CellType==j,])[1])
        temp_pd = pd[complete.cases(pd),]
        temp_beta = beta[,colnames(beta)%in%rownames(temp_pd)]
        design = model.matrix(~temp_pd$status)
        fit = lmFit(beta2m(temp_beta), design)
        fit = eBayes(fit)
        table_limma = topTable(fit, coef=2, n=dim(fit)[1])
        # table_limma = topTable(fit, coef=2, n=dim(fit)[1], adjust.method="BH")
        table_limma_complete = table_limma[is.finite(rowSums(table_limma)),]
        sig_limma = table_limma_complete[table_limma_complete$adj.P.Val<p_val,]
        hypo_limma = sig_limma[order(sig_limma[, "logFC"], decreasing = FALSE), ]
        hyper_limma = sig_limma[order(sig_limma[, "logFC"], decreasing = TRUE), ]
        # probes = c(rownames(hyper_limma)[seq_len(floor(n_prob/2))],
        #            rownames(hypo_limma)[seq_len(floor(n_prob/2))])
        probes = c(rownames(sig_limma))
        sig_limma_list = append(sig_limma_list,list(probes))
      }
    }
    newnewSetofMarkers = append(newnewSetofMarkers,list(Reduce(union, sig_limma_list)))
  }
  newnewSetofMarkers = unique(unlist(newnewSetofMarkers))
  newnewSetofMarkers = newnewSetofMarkers[complete.cases(newnewSetofMarkers)]
  return((save(newnewSetofMarkers,file="newnewSetofMarkers.RData")))
}

# find the markers:
beta <- comb_beta
pd <- as.data.frame(comb_pd)
celltypes <- c('Bcell','CD4T','CD8T','Gran','Mono','NK','nRBC')
p_val <- 1e-8
# p_val <- 0.05

deconv_pair(beta, pd, celltypes, p_val)

# get the new set of markers' matrix:
load("/home/liuwent/04b-cell_type_deconvolution/newnewSetofMarkers.RData")
length(newnewSetofMarkers)
# 381,923 or 151,794

newnewSetofMarkers_pd <- comb_pd
dim(newnewSetofMarkers_pd)
# 87   2
# save(newnewSetofMarkers_pd, file="newnewSetofMarkers_pd.RData")
# load("/home/liuwent/04b-cell_type_deconvolution/newnewSetofMarkers_pd.RData")

newnewSetofMarkers_beta <- comb_beta[which(rownames(comb_beta)%in%newnewSetofMarkers),]
dim(newnewSetofMarkers_beta)
# 381,923 or 151,794   87
# save(newnewSetofMarkers_beta, file="newnewSetofMarkers_beta.RData")
# load("/home/liuwent/04b-cell_type_deconvolution/newnewSetofMarkers_beta.RData")

## PCA check:
beta = t(newnewSetofMarkers_beta)
beta_stand = sweep(beta,2,colMeans(beta),"-")
SVD <- svd(beta_stand)
svd_plot = cbind(as.data.frame(rownames(newnewSetofMarkers_pd)),
                 as.data.frame(SVD$u[, 1]), 
                 as.data.frame(SVD$u[, 2]), 
                 as.data.frame(newnewSetofMarkers_pd$CellType))
colnames(svd_plot) = c("ID","PC1","PC2","CellType")
col.pat = brewer.pal(n = length(unique(newnewSetofMarkers_pd$CellType)), name = "Set2")

pdf("/home/liuwent/04b-cell_type_deconvolution/newnewSetofMarkers_PCA.pdf")
ggplot(svd_plot, aes(x=PC1, y=PC2, color=CellType)) +
  scale_color_manual(values=col.pat)+
  geom_point()+
  theme_classic(base_size = 18)
dev.off()


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
# save(newnewSetofMarkers_final, file = "newnewSetofMarkers_final.RData")
# load("/home/liuwent/04b-cell_type_deconvolution/newnewSetofMarkers_final.RData")



