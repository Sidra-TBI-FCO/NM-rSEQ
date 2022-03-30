# Setup environment
rm(list=ls())

setwd("./")

# Load libraries 
library(stringr)
library(ggplot2)
library(openxlsx)

# Set parameters
cutoff = "0.01"  # 0.01 # 0.05

# Load data
IPA = read.xlsx(paste0("./raw_data_intersect_up_down_EDAseq_",cutoff,".xlsx"))

# Number of downregulated genes 
IPA$Downregulated = gsub(".*\\(", "", IPA$Downregulated)
IPA$Downregulated = gsub("%)", "", IPA$Downregulated)

# Number of upregulated genes 
IPA$Upregulated = gsub(".*\\(", "", IPA$Upregulated)
IPA$Upregulated = gsub("%)", "", IPA$Upregulated)

# Upregulated pathways
IPA.up = data.frame(pathway = IPA$Ingenuity.Canonical.Pathways)
IPA.up$direction = "Up"
IPA.up$pct = IPA$Upregulated[match(IPA.up$pathway, IPA$Ingenuity.Canonical.Pathways)]
IPA.up$logp = IPA$`-log(p-value)`[match(IPA.up$pathway, IPA$Ingenuity.Canonical.Pathways)]
IPA.up$zscore = IPA$`z-score`[match(IPA.up$pathway, IPA$Ingenuity.Canonical.Pathways)]

# Downregulated pathways
IPA.down = data.frame(pathway = IPA$Ingenuity.Canonical.Pathways)
IPA.down$direction = "down"
IPA.down$pct = IPA$Downregulated[match(IPA.down$pathway, IPA$Ingenuity.Canonical.Pathways)]
IPA.down$logp = IPA$`-log(p-value)`[match(IPA.down$pathway, IPA$Ingenuity.Canonical.Pathways)]
IPA.down$zscore = IPA$`z-score`[match(IPA.down$pathway, IPA$Ingenuity.Canonical.Pathways)]

# Combine upregulated and downregulated pathways 
final.IPA = rbind(IPA.up, IPA.down)
final.IPA$group = "Intersect Genes"

# Order based on log pvalue
final.IPA = final.IPA[order(final.IPA$logp, decreasing = F),]

# Set pathways as factor 
final.IPA$pathway = factor(final.IPA$pathway, levels = unique(final.IPA$pathway))
pw_label = levels(final.IPA$pathway)

# Set pct as numeric
final.IPA$pct = as.numeric(final.IPA$pct)

# plot 
svg(filename = paste0("./intersect_up_down_genes_EDAseq_",cutoff,".svg"), height = 5,width = 7) 

ggplot(data = final.IPA,aes(x=pathway,group=group,fill=direction)) +
  geom_bar(aes(y=pct),position="stack", stat="identity") +
  scale_fill_manual(values = c("green","red")) +
  geom_line(aes(y=logp*50/10),color="orange",size=0.5) + 
  geom_point(aes(y=logp*50/10),color="orange",size=1) +  
  scale_y_continuous(name = "Enrichment(%)",
                     sec.axis = sec_axis(~ . * 10 / 50 , name = "-logP10"), 
                     limits = c(-3, 50)) +
  geom_hline(yintercept=-log(0.05,10)*50/10, linetype="dashed", color = "blue") +  
  geom_hline(yintercept=-log(0.005,10)*50/10, linetype="dashed", color = "blue", size = 0.5) + 
  geom_point(aes(y=-3,color=zscore),size=3) +
  scale_color_gradient2(low = "blue", mid = "white", high = "#FF8C00", na.value = "#C0C0C0") +
  scale_x_discrete(labels=pw_label) +
  theme_bw() +
  facet_wrap( ~ group,nrow = 1,ncol = 2) +
  coord_flip()

dev.off()

