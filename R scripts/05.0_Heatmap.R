# Selected rows labeled 

# All DEG in all conditions 

# Setup environment
rm(list=ls())

setwd("./")

# load libraries 
library(ComplexHeatmap)
library(stringr)
library(openxlsx)
library(circlize)

# Load data 
load("./EDAseq_normalized_selected_sample.Rdata")
Annotation = read.csv("./Annotation_file.csv")
DEG = read.csv("./all_DEG_FDR_0.01.csv") # csv file was ordered based on gene names from A to z prior loading 

data_nor = data_log2

# subset counts
data_nor = data_nor[which(rownames(data_nor) %in% DEG$x),]

# Subset sample info
Annotation = Annotation[which(Annotation$Condition %in% c("CTRL", "LPS", "ConA", "MX1", "MX3", "MX5")),]

# Subset counts matrix
data_nor = data_nor[,which(colnames(data_nor) %in% Annotation$sampleID_original)]

colnames(data_nor) == Annotation$sampleID_original

Annotation$Condition = factor(Annotation$Condition, levels = c("CTRL", "LPS", "ConA", "MX1", "MX3", "MX5"))

table(Annotation$MATERIAL)
Annotation$MATERIAL = factor(Annotation$MATERIAL, levels = c("CTRL", "LPS", "ConA", "Ta4C3", "Mo2Ti2C3", "Nb4C3"))

Expression.matrix = data_nor
#Expression.matrix = log(Expression.matrix+1,2)

class(Expression.matrix)
Expression.matrix = as.matrix(Expression.matrix)
# z-score Expression.matrix
Expression.matrix.z = Expression.matrix

for(j in 1: nrow(Expression.matrix.z))  {
  Expression.matrix.z[j,] = (Expression.matrix[j,]-mean(Expression.matrix[j,]))/sd(Expression.matrix[j,]) # z-score the enrichment matrix
}

colnames(Expression.matrix.z) == Annotation$sampleID_original

ha = HeatmapAnnotation(df = data.frame(Material = Annotation$MATERIAL),
                       col = list(Material = c("CTRL" = "#000000", "LPS" = "#808080", "ConA" = "#20B2AA", 
                                               "Ta4C3" = "#709ABB", "Mo2Ti2C3" = "#8CC8F7", "Nb4C3" = "#A8F9E2")),
                       show_legend = TRUE)

# setting z score scale values and the colors (min = minumum, max = maximum)
colors = colorRamp2(c(min(Expression.matrix.z), 0, max(Expression.matrix.z)), c("black", "blue", "yellow"))

r.n = as.data.frame(rownames(Expression.matrix.z))
labels = rownames(Expression.matrix.z)[c(1879,1896,1897,1919)]

ha.r = rowAnnotation(Expression.matrix.z = anno_mark(at = c(1879,1896,1897,1919), labels = rownames(Expression.matrix.z)[c(1879,1896,1897,1919)]))

svg(filename = paste0("./limma_FDR_0.01_all_DEG_new.svg"), height = 5,width = 6)
genes.heatmap = Heatmap(Expression.matrix.z,
                        name = "Expression\n z score",
                        cluster_rows = TRUE,
                        cluster_columns = T,
                        show_heatmap_legend = TRUE,
                        top_annotation = ha,
                        row_title_gp = gpar(fontsize = 2),
                        column_names_gp = gpar(fontsize = 4),
                        row_names_gp = gpar(fontsize = 2),
                        show_column_names = F,
                        show_row_names = F,
                        right_annotation = ha.r,
                        col = colors, 
                        row_names_max_width = unit(8, "in"))

draw(genes.heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

#################################################################################################################################################
#################################################################################################################################################
# DEG intersect in MX

# Setup environment
rm(list=ls())

setwd("./")

# load libraries 
library(ComplexHeatmap)
library(stringr)
library(openxlsx)
library(circlize)

# Load data 
load("./EDAseq_normalized_selected_sample.Rdata")
Annotation = read.csv("./Annotation_file.csv")
DEG = read.csv("./intersect_up_down_FDR_0.01.csv")

data_nor = data_log2

# subset counts
data_nor = data_nor[which(rownames(data_nor) %in% DEG$x),]

# Subset sample info, expression matrix to selected samples
Annotation = Annotation[which(Annotation$Condition %in% c("CTRL", "MX1", "MX3", "MX5")),]

# Subset counts matrix
data_nor = data_nor[,which(colnames(data_nor) %in% Annotation$sampleID_original)]

colnames(data_nor) == Annotation$sampleID_original

Annotation$Condition = factor(Annotation$Condition, levels = c("CTRL", "MX1", "MX3", "MX5"))

table(Annotation$MATERIAL)
Annotation$MATERIAL = factor(Annotation$MATERIAL, levels = c("CTRL", "Ta4C3", "Mo2Ti2C3", "Nb4C3"))

Expression.matrix = data_nor

class(Expression.matrix)
Expression.matrix = as.matrix(Expression.matrix)

# z-score Expression.matrix
Expression.matrix.z = Expression.matrix

for(j in 1: nrow(Expression.matrix.z))  {
  Expression.matrix.z[j,] = (Expression.matrix[j,]-mean(Expression.matrix[j,]))/sd(Expression.matrix[j,]) # z-score the enrichment matrix
}

colnames(Expression.matrix.z) == Annotation$sampleID_original

ha = HeatmapAnnotation(df = data.frame(Material = Annotation$MATERIAL),
                       col = list(Material = c("CTRL" = "#000000", "Ta4C3" = "#709ABB", "Mo2Ti2C3" = "#8CC8F7", "Nb4C3" = "#A8F9E2")),
                       show_legend = TRUE)

# setting z score scale values and the colors (min = minumum, max = maximum)
colors = colorRamp2(c(min(Expression.matrix.z), 0, max(Expression.matrix.z)), c("black", "blue", "yellow"))


svg(filename = paste0("./imma_intersect_up_down_MXs_FDR_0.01.svg"), height = 5,width = 6)
genes.heatmap = Heatmap(Expression.matrix.z,
                        name = "Expression\n z score",
                        cluster_rows = TRUE,
                        cluster_columns = F,
                        show_heatmap_legend = TRUE,
                        top_annotation = ha,
                        row_title_gp = gpar(fontsize = 3),
                        column_names_gp = gpar(fontsize = 3),
                        row_names_gp = gpar(fontsize = 8),
                        show_column_names = F,
                        col = colors,
                        row_names_max_width = unit(7, "in"))

draw(genes.heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

