# Setup environment

rm(list=ls())

setwd("./")

# Set Parameters 
cutoff = "0.01"  # 0.01  0.05 

# Load data 
load("./EDAseq_normalized_selected_sample.Rdata")
Annotation = read.csv("./Annotation_file.csv")
DEG = read.csv("./MX1_CTRL_all_DEG_EDAseq.csv")
intersect = read.csv(paste0("./intersect_up_down_FDR_",cutoff,".csv"))

FC = DEG[which(DEG$X %in% intersect$genes),]
write.csv(FC, file = paste0("./FC_for_intersct_up_down_",cutoff,".csv"))
