rm(list=ls())
setwd("~/KayleeStuff/Smith_Lab/Data/RNA_Seq/Mat_Liv")
library("FactoMineR")
library("ggplot2")
library("factoextra")
library("tibble")
library("ggsci")

#Import table of differentially expressed genes and remove extraneous columns
DEG <- read.table("DEG_Results_MDvAlc_MatLiv_LFCShrink_VSTParam_RemovedNA_ExcelProof.txt", header = TRUE, sep = "", row.names = "GeneID")
DEG_reduced <- DEG[,7:22]
DEG_reduced_flip <- t(DEG_reduced) #reverse rows and columns
group <- factor(c(rep("Control",8), rep("Alcohol",8))) #create group name list
group <- data.frame(group) #turn group name list into dataset
DEG_reduced_flip_labeled <- cbind.data.frame(DEG_reduced_flip, group) #add in group name column
DEG_reduced_flip_labeled <- DEG_reduced_flip_labeled[,c(14789, 1:14788)] #rearrange the columns to put group names 1st

#Create PCA but do not graph it. 
results.pca <- PCA(DEG_reduced_flip_labeled[,-1], scale.unit = FALSE, ncp = 5, graph = FALSE) #ncp = number of principal components
print(results.pca)

#Create scree plot
png("MatLiv_DEG_ScreePlot_highRes.png", units="in", width=7, height=7, res=600) 
fviz_eig(results.pca, addlabels = TRUE, ylim = c(0,90))
dev.off()

#See which variables contribute the most to each of the principal components
variables <- get_pca_var(results.pca)
variables$coord
contributions <- head(variables$contrib, 20)

#Look at individual contributions and graph PCA and export the graph
png("MatLiv_DEG_PCA_basic_highRes.png", units="in", width=9, height=7, res=600) 
ind.pca <- fviz_pca_ind(results.pca, col.ind = DEG_reduced_flip_labeled$group, palette = "aaas", addEllipses = TRUE, legend.title = "Group", mean.point = FALSE, repel = TRUE)
ggpubr::ggpar(ind.pca, title = "Principal Component Analysis", subtitle = "Maternal Liver") #, xlab = "PC1 (61.3%)", ylab = "PC2 (24.8%)"
dev.off()

#alt PCA
png("MatLiv_DEG_PCA_EllipsesConfidence_highRes.png", units="in", width=9, height=7, res=600) 
ind.pca <- fviz_pca_ind(results.pca, col.ind = DEG_reduced_flip_labeled$group, palette = "aaas", addEllipses = TRUE, ellipse.type = "confidence", slegend.title = "Group", mean.point = FALSE, repel = TRUE)
ggpubr::ggpar(ind.pca, title = "Principal Component Analysis", subtitle = "Maternal Liver") #, xlab = "PC1 (61.3%)", ylab = "PC2 (24.8%)"
dev.off()

#PCA without MD4
DEG_reduced_flip_noMD4 <- DEG_reduced_flip[c(1:3,5:16),]
group2 <- factor(c(rep("Control",7), rep("Alcohol",8))) #create group name list
group2 <- data.frame(group2) #turn group name list into dataset
DEG_reduced_flip_labeled_noMD4 <- cbind.data.frame(DEG_reduced_flip_noMD4, group2) #add in group name column
DEG_reduced_flip_labeled_noMD4 <- DEG_reduced_flip_labeled_noMD4[,c(14789, 1:14788)] #rearrange the columns to put group names 1st

results.pca.noMD4 <- PCA(DEG_reduced_flip_labeled_noMD4[,-1], scale.unit = FALSE, ncp = 5, graph = FALSE) #ncp = number of principal components
png("MatLiv_DEG_PCA_NoMD4_highRes.png", units="in", width=9, height=7, res=600) 
ind.pca2 <- fviz_pca_ind(results.pca.noMD4, col.ind = DEG_reduced_flip_labeled_noMD4$group, palette = "aaas", addEllipses = TRUE, legend.title = "Group", mean.point = FALSE, repel = TRUE)
ggpubr::ggpar(ind.pca2, title = "Principal Component Analysis", subtitle = "Maternal Liver") #, xlab = "PC1 (61.3%)", ylab = "PC2 (24.8%)"
dev.off()

png("MatLiv_DEG_PCA_NoMD4_EllipsesConfidence_highRes.png", units="in", width=9, height=7, res=600) 
ind.pca2 <- fviz_pca_ind(results.pca.noMD4, col.ind = DEG_reduced_flip_labeled_noMD4$group, palette = "aaas", addEllipses = TRUE, ellipse.type = "confidence", legend.title = "Group", mean.point = FALSE, repel = TRUE)
ggpubr::ggpar(ind.pca2, title = "Principal Component Analysis", subtitle = "Maternal Liver") #, xlab = "PC1 (61.3%)", ylab = "PC2 (24.8%)"
dev.off()

#HCPC analysis
res.hcpc <- HCPC(results.pca, graph = FALSE)
png("MatLiv_DEG_HCPC_basic_highRes.png", units="in", width=9, height=10, res=600) 
fviz_dend(res.hcpc, cex = 0.7, palette = "aaas", rect = TRUE, rect_fill = TRUE, rect_border = "aaas", labels_track_height = 0.8)
dev.off()
