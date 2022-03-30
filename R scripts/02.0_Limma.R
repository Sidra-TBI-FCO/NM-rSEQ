# Setup environment
rm(list=ls())

setwd("./")

# load libraries 
library(limma)
library(ggplot2)
library(ggpubr)
library(stringr)

# load data
load("./EDAseq_normalized_selected_sample.Rdata")
Annotation = read.csv("./Annotation_file.csv")

nor.matrix = as.matrix(data_log2)

Annotation$sampleID_original == colnames(nor.matrix)

table(Annotation$MATERIAL)

# Set Parameters 
cond1 = "CTRL"
cond2 = "LPS" # LPS ConA MX1 MX3 MX5  

annotation = Annotation
counts.matrix = nor.matrix

annotation = annotation[which(Annotation$Condition %in% c(cond1, cond2)),]
counts.matrix = counts.matrix[,which(colnames(counts.matrix) %in% annotation$sampleID_original)]

colnames(counts.matrix) == annotation$sampleID_original

# Limma 
dim(counts.matrix)
counts.matrix = counts.matrix[which(rowSums(counts.matrix) > 0),] 

dim(counts.matrix)

counts.matrix = as.data.frame(t(counts.matrix))
counts.matrix$group = NA
counts.matrix$group = annotation$Condition[match(rownames(counts.matrix), annotation$sampleID_original)]
counts.matrix$group = factor(counts.matrix$group, levels = c(paste0(cond1), paste0(cond2)))

design=model.matrix(~ counts.matrix$group)
counts.matrix$group = NULL

counts.matrix = as.matrix(counts.matrix)
mode(counts.matrix) = "numeric"

fit = lmFit(t(counts.matrix) , design=design)
fit = eBayes(fit)

diff.stats = topTable(fit, coef=2,adjust.method = "fdr", number = ncol(counts.matrix))
write.csv(diff.stats, file = paste0("./",cond2,"_", cond1,"_all_DEG_EDAseq.csv"))

# 0.01 
DEG1 = diff.stats[which(diff.stats$adj.P.Val < 0.01),]
length(rownames(DEG1))
write.csv(DEG1, file = paste0("./",cond2,"_", cond1,"_FDR_0.01_DEG_EDAseq.csv"))

# 0.05
DEG2 = diff.stats[which(diff.stats$adj.P.Val < 0.05),]
length(rownames(DEG2))
write.csv(DEG2, file = paste0("./",cond2,"_", cond1,"_FDR_0.05_DEG_EDAseq.csv"))

# 0.1
DEG3 = diff.stats[which(diff.stats$adj.P.Val < 0.1),]
length(rownames(DEG3))
