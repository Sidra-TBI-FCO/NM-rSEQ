# Setup environment
rm(list=ls())

setwd("./")

# load libraries 
library(ggplot2)
library(stringr)
library(openxlsx)

genes = read.xlsx("./numbers_DEG_EADseq.xlsx")


FDR1 = genes[,c(8,9)]
colnames(FDR1)[2] = "Nubmer of Genes"
FDR1$p.adj = "FDR < 0.01"

FDR2 = genes[,c(8,10)]
colnames(FDR2)[2] = "Nubmer of Genes"
FDR2$p.adj = "FDR < 0.05"

FDR3 = genes[,c(8,11)]
colnames(FDR3)[2] = "Nubmer of Genes"
FDR3$p.adj = "FDR < 0.1"

merge = rbind(FDR1, FDR2)

FDR = rbind(merge, FDR3)

FDR$MATERIAL = factor(FDR$MATERIAL, levels = c("ConA", "LPS", "Nb4C3", "Mo2Ti2C3", "Ta4C3"))
FDR$p.adj = factor(FDR$p.adj, levels = c("FDR < 0.01", "FDR < 0.05", "FDR < 0.1"))


svg(filename = paste0("./numbers_DEG_MATERIAL_EDAseq.svg"), height = 4,width = 7)

ggplot(FDR, aes(fill=p.adj, y=`Nubmer of Genes`, x=MATERIAL)) + 
  geom_bar(position="dodge", stat="identity") +
  ylim(0,8000) +
  geom_text(aes(label = signif(`Nubmer of Genes`, digits = 3)), nudge_y = 0.1, cex = 2.5) +
  scale_fill_manual(values=c("#FFDAB9", "#AFEEEE", "#C0C0C0")) + theme(text = element_text(size = 8), axis.text = element_text(size = 8)) + 
theme_classic() 

dev.off()
