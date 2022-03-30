# Setup environment
rm(list=ls())

setwd("./")

# Load libraries 
library(EDASeq)
library(base64enc)
library(preprocessCore)

# Load required files
load("./geneInfo.July2017.RData")
counts = "./DB_SDR400025_201202_LEX8_B1_Mexenes_raw_counts.txt"
sub.reads = readLines(counts, n=1) 
raw.counts = read.csv(counts, sep = "\t",  stringsAsFactors = F)

samples.data = raw.counts

samples.data = aggregate(samples.data,FUN="mean", by=list(samples.data$symbol))

#Edit samples.data format
#Rownames 
rownames(samples.data) = samples.data$Group.1

samples.data$ID = NULL
samples.data$symbol = NULL
samples.data$Group.1 = NULL

#convert samples.data into matrix 
samples.matrix = as.matrix(samples.data)

# Set geneInfo as data.freme 
geneInfo = as.data.frame(geneInfo)

#Extract genes from geneInfo that matches genes in samples.matrix
common.genes = unique(rownames(samples.matrix)[which(rownames(samples.matrix) %in% rownames(geneInfo))])

# Extract the common genes for geneInfo and samples.matrix
geneInfo = geneInfo[common.genes,]
samples.filtered = samples.matrix[common.genes,]

mode(samples.filtered) = "numeric"
dim(samples.filtered) 

# Data Normalization using EDASeq
samples.exp.norm = newSeqExpressionSet(samples.filtered, featureData = geneInfo)

# removes effect related to between lane distributional differences, as sequencing depth
samples.exp.norm = betweenLaneNormalization(samples.exp.norm, which = "upper", offset = T)

#Take= log (unnormalized + .1) + offst(normalized)
samples.norm.log = log(samples.filtered +.1) + offst(samples.exp.norm)
samples.norm.log = floor(exp(samples.norm.log) - .1)  

#Quantile Normalization
samples.quantiles.norm = normalize.quantiles(samples.norm.log)
samples.quantiles.norm = floor(samples.quantiles.norm)

rownames(samples.quantiles.norm) = rownames(samples.norm.log)
colnames(samples.quantiles.norm) = colnames(samples.norm.log)

#Log transformation (final normalized-transformed file)
data_log2 = log(samples.quantiles.norm+1,2) #log base 2
dim(data_log2)


save(data_log2, file="./EDAseq_normalized_selected_sample.Rdata")


