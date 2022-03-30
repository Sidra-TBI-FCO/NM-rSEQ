# intersect with all materials 

# Setup environment
rm(list=ls())

setwd("./")

# load libraries 
library(VennDiagram)
library(venn)
library(ggplot2)

# Set parameters 
cutoff = "0.01"  # 0.01  # 0.05

# Load data
LPS = read.csv(paste0("./LPS_CTRL_FDR_",cutoff,"_DEG_EDAseq.csv"))
ConA = read.csv(paste0("./ConA_CTRL_FDR_",cutoff,"_DEG_EDAseq.csv"))
MX1 = read.csv(paste0("./MX1_CTRL_FDR_",cutoff,"_DEG_EDAseq.csv"))
MX3 = read.csv(paste0("./MX3_CTRL_FDR_",cutoff,"_DEG_EDAseq.csv"))
MX5 = read.csv(paste0("./MX5_CTRL_FDR_",cutoff,"_DEG_EDAseq.csv"))


up.LPS = LPS[which(LPS$logFC > 0),]
up.ConA = ConA[which(ConA$logFC > 0),]
up.MX1 = MX1[which(MX1$logFC > 0),]
up.MX3 = MX3[which(MX3$logFC > 0),]
up.MX5 = MX5[which(MX5$logFC > 0),]

down.LPS = LPS[which(LPS$logFC < 0),]
down.ConA = ConA[which(ConA$logFC < 0),]
down.MX1 = MX1[which(MX1$logFC < 0),]
down.MX3 = MX3[which(MX3$logFC < 0),]
down.MX5 = MX5[which(MX5$logFC < 0),]

# Intersect MXs
# Downregulated
venn.plot = 
  list(Ta4C3 = down.MX1$X, Mo2Ti2C3 = down.MX3$X, Nb4C3 = down.MX5$X)

svg("./limma_down_FDR_0.01_MXs.svg", width = 5, height = 5)

venn.result =
  venn(venn.plot, ilabels = TRUE, 
       zcolor = c("#709ABB", "#8CC8F7", "#A8F9E2"), size = 100, cexil = 2, 
       cexsn = 10, ilcs = 2, sncs = 2);

dev.off()


# Upregulated
venn.plot = 
  list(Ta4C3 = up.MX1$X, Mo2Ti2C3 = up.MX3$X, Nb4C3 = up.MX5$X)

svg("./limma_up_FDR_0.01_MXs.svg", width = 5, height = 5)

venn.result =
  venn(venn.plot, ilabels = TRUE, 
       zcolor = c("#709ABB", "#8CC8F7", "#A8F9E2"), size = 100, cexil = 2, 
       cexsn = 10, ilcs = 2, sncs = 2);


dev.off()


#############################################

# Save intersect genes 

common.up = as.data.frame(Reduce(intersect,list(up.MX1$X, up.MX3$X, up.MX5$X)))
common.down = as.data.frame(Reduce(intersect,list(down.MX1$X, down.MX3$X, down.MX5$X)))

colnames(common.down)[1] = "genes"
colnames(common.up)[1] = "genes"

intersect.genes = rbind(common.down, common.up)

write.csv(intersect.genes$genes, file= paste0("./intersect_up_down_FDR_",cutoff,".csv"))


# All DEG

all.up = rbind(up.ConA, up.LPS, up.MX1, up.MX3, up.MX5)
all.down = rbind(down.ConA, down.LPS, down.MX1, down.MX3, down.MX5)

colnames(all.up)[1] = "genes"
colnames(all.down)[1] = "genes"

all.DEG = rbind(all.up, all.down)

write.csv(all.DEG$genes, file = paste0("./all_DEG_FDR_",cutoff,".csv"))
