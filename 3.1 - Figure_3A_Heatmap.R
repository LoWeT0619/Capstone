# Library:
library(dplyr)
library(limma)
library(ggplot2)
library(tidyverse)
library(gplots)

# Data:
load("/home/liuwent/04-Full_Model/pd.RData")
load("/home/liuwent/04b-cell_type_deconvolution/estF.RData")

estF <- as.data.frame(estF)
rownames(estF) <- paste(pd$Sample_Name, pd$Sample_Group, sep = "-")

# Heatmap.2:
pdf("/home/liuwent/04b-cell_type_deconvolution/Fig_3A_heatmap.pdf")
col_fun <- colorRampPalette(c("blue", "green", "yellow", "red"))

heatmap.2(as.matrix(estF), 
          trace = "none", 
          dendrogram = "both", 
          col = col_fun(100),
          scale = "none")

custom_colors <- col_fun(100)
custom_breaks <- c(-Inf, 0, 0.05, 0.1, 0.3, 0.5, Inf)
custom_values <- seq(0, 1, length.out = length(custom_breaks))

heatmap_scale <- scale_color_gradientn(colors = custom_colors, 
                                       values = custom_values, 
                                       breaks = custom_breaks)

heatmap.2(as.matrix(estF), 
          trace = "none", 
          dendrogram = "both", 
          col = col_fun(100),
          scale = "none") + 
  heatmap_scale
dev.off()

# Heatmap - ggplot2:
hm_data <- estF %>% 
  rownames_to_column() %>%
  gather(colname, value, -rowname) %>% 
  setNames(c("Samples", "Cell_Type_Proportion", "Value"))

pdf("/home/liuwent/04b-cell_type_deconvolution/Fig_3A_heatmap.pdf")
ggplot(data = hm_data, aes(x=Cell_Type_Proportion, y=Samples, fill=Value)) +
  geom_tile() + 
  scale_fill_gradient2(low="#0000FF", mid="#FFFFFF", high="#FF0000")
dev.off()

