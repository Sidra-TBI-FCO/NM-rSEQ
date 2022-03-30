# Setup environment
rm(list=ls())

setwd("./")

# load libraries 
library(scatterplot3d)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(ggrepel)
library(openxlsx)

# Load data
load("./EDAseq_normalized_selected_sample.Rdata")
Annotation = read.csv("./Annotation_file.csv")

Annotation$sampleID_original == colnames(data_log2)

table(Annotation$MATERIAL)
Annotation$MATERIAL = factor(Annotation$MATERIAL, levels = c("CTRL", "LPS", "ConA", "Ta4C3", "Mo2Ti2C3", "Nb4C3"))

colors2 = c("CTRL" = "#000000", "LPS" = "#808080", "ConA" = "#20B2AA", "Ta4C3" = "#709ABB", "Mo2Ti2C3" = "#8CC8F7", "Nb4C3" = "#A8F9E2")

QN.log2 = t(data_log2)

quant.Log2.transformed = QN.log2
dim(quant.Log2.transformed)

quant.zero.genes.removed = quant.Log2.transformed[,apply(quant.Log2.transformed,2,var, na.rm = T) !=0]

dim(quant.zero.genes.removed)

# Principal Component Analysis
PCA = prcomp(quant.zero.genes.removed, center = T, scale = T)
percent.var = round(100*(PCA$sdev^2 / sum(PCA$sdev^2)),1)

# Extract the pca scores 
pca.scores = data.frame(PCA$x[,1:3])
pca.scores = pca.scores[match(Annotation$sampleID_original, rownames(pca.scores)),]

# Add "variable" column 
colored.by = "MATERIAL"

pca.scores$variable = Annotation[,colored.by]
pca.scores$variable = as.factor(pca.scores$variable)

# Add each color to the corresponding type for each sample
condition.color = colors2[as.numeric(pca.scores$variable)]

# 2D PCA
svg(filename = paste0("./PCA_materials_EDAseq_new.svg"), height = 4,width = 6)

ggplot(pca.scores, aes(x=PC1, y=PC2, color = variable)) +
  geom_point(size=4, aes(col=variable))+
  scale_color_manual(name="Materials", values = colors2) +
  labs(x=paste0("PC1 (",percent.var[1],"%)"),y=paste0("PC2 (",percent.var[2],"%)"))+
  stat_ellipse()+  
  theme_light(base_size=20) +
  theme(axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.title = element_text(size = 15))
dev.off()
