rm(list=ls())
setwd("~/KayleeStuff/Smith_Lab/Data/RNA_Seq/Fet_Liv")
library("DESeq2")
library("ggplot2")
library("ggpubr")
library("DEP")
library("hexbin")
library("devtools")
library("apeglm")
library("pheatmap")
library("RColorBrewer")
library("SummarizedExperiment")
library("genefilter")
library("BiocManager")

#Input file, select out correct file columns, rename file columns, turn table into matrix, create treatment table for colData 
genes <- read.table("FetMLiver_MdvsEtOH.txt", header = TRUE, sep = "", row.names = "Geneid")
genes.reduced <- genes[,c(6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)]
names(genes.reduced) <- c("MD-E1", "MD-E2", "MD-E3", "MD-E4", "MD-E5", "MD-E6", "MD-E7", "MD-E8", 
                          "Alc-F1", "Alc-F2", "Alc-F4", "Alc-F5", "Alc-F6", "Alc-F7", "Alc-F8", "Alc-F9")
countdata <- as.matrix(genes.reduced) #Table must be converted to a matrix for DeSeq to perform calculations
treatment <- factor(c(rep("Control",8), rep("Alcohol",8)))
coldata <- data.frame(row.names = colnames(countdata), treatment)

#Ensure that the rownames of the matrix are in the same order as in the other file. DESeq will not check. 
#If correct, running this code will output "TRUE". 
all (rownames(colData) == colnames(genes.reduced))

#Make a DESeq data set from the count matrix and treatment information
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~treatment) #This step gives an error, but it is corrected in the next line of code
dds$treatment <- factor(dds$treatment, levels = c("Control","Alcohol")) #sets the comparison of alcohol vs. control for DeSeq2 step

#If you want to filter out the reads that have low copy number, then you can run this code
keep <- rowSums(counts(dds) >= 10) >=5
table(keep)
dds <- dds[keep,]

##Perform both the rlog and the vst transformation and then compare the meanSdPlot and PCA for them. 
#Make sure to export the resulting datasets from these steps for further PCA and other analyses. 
#For decreasing gene-wide dispersion use parametric if it isn't influenced by the data as much as local or mean. Find the trend that is least affected by the dispersion

#rlogTransformation
#rlog_parametric <- rlogTransformation(dds, blind = FALSE, fitType = "parametric") 
#rlog_local <- rlogTransformation(dds, blind = FALSE, fitType = "local")
#rlog_mean <- rlogTransformation(dds, blind = FALSE, fitType = "mean")

#head(assay(rlog_mean), 10)
#colData(rlog_mean)
#meanSdPlot <- meanSdPlot(rlog_mean)
#plotPCA(rlog_mean, intgroup = c("treatment"))

#vstTransformation
vst_parametric <- varianceStabilizingTransformation(dds, blind = FALSE, fitType = "parametric")
#vst_local <- varianceStabilizingTransformation(dds, blind = FALSE, fitType = "local")
#vst_mean <- varianceStabilizingTransformation(dds, blind = FALSE, fitType = "mean")

head(assay(vst_parametric), 10)
colData(vst_parametric)
meanSdPlot2 <- meanSdPlot(vst_parametric)
plotPCA(vst_parametric, intgroup = c("treatment"))

#P-value calculation and DeSeq2 step
dds <- DESeq(dds)

#Plot gene dispersions
plotDispEsts(dds)

#Results
#To get the correct coefficient, run "resultsNames(dds)" to figure out what to put after "coef="". 
#In "contrast" it should be: the name of the factor in the design formula, then the control, then the treated
results_unshrunken <- results(dds)
results_lfc <- lfcShrink(dds, coef="treatment_Alcohol_vs_Control",type="apeglm")
plotMA(results_unshrunken, ylim=c(-2,2))
plotMA(results_lfc, ylim=c(-2,2))
#order results (res) by pvalue
results_Ordered <- results_lfc[order(results_lfc$padj), ]

#summarize the results
summary(results_Ordered)
sum(results_unshrunken$padj < 0.1, na.rm=TRUE)

#merge results data frame with counts data frame to create a table that contains both counts and statistical results
results_lfc_combined <- merge(as.data.frame(results_lfc), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(results_lfc_combined)[1] <- "GeneID"
results_lfc_combined_rmNA <- na.omit(results_lfc_combined) #remove genes with NA values
results_unshrunken_combined <- merge(as.data.frame(results_unshrunken), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(results_unshrunken_combined)[1] <- "GeneID"
results_unshrunken_combined_rmNA <- na.omit(results_unshrunken_combined) #remove genes with NA values

#calculate the number of DE (true) and non-DE (false) genes
deCount <- table(results_lfc$padj<0.05)

#write to CSV
write.table(results_lfc_combined, file="DEG_Results_MDvAlc_FetLiv_LFCShrink_VSTParam.txt", sep = "\t")
write.table(results_lfc_combined_rmNA, file="DEG_Results_MDvAlc_FetLiv_LFCShrink_VSTParam_RemovedNA.txt", sep = "\t")
write.table(results_unshrunken_combined, file="DEG_Results_MDvAlc_FetLiv_NoShrink_VSTParam.txt", sep = "\t")
write.table(results_unshrunken_combined_rmNA, file="DEG_Results_MDvAlc_FetLiv_NoShrink_VSTParam_RemovedNA.txt", sep = "\t")

#Plot the counts for a single normalized gene
plotCounts(dds, "Cebpa", intgroup = "treatment", normalized = TRUE, main = "Cebpa")

#Histogram of p-values
hist(results_lfc$pvalue, breaks=50, col="grey")

#Pval adjustment using BH and Bonferonni methods as a check of methods
adjusted_pvals_BH <- p.adjust(results_lfc$pvalue, method="BH")
adjusted_pvals_bon <- p.adjust(results_lfc$pvalue, method="bonferroni")
hist(adjusted_pvals_BH)
hist(adjusted_pvals_bon)
BHtotal <- sum(adjusted_pvals_BH < 0.05, na.rm=TRUE) #removes NA values
Bontotal <- sum(adjusted_pvals_bon < 0.1, na.rm=TRUE)

#Heatmap for sample distance
sampleDists <- dist(t(assay(vst_parametric)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette((rev(brewer.pal(9, "Blues"))))(255)
#Heatmap with the rows and columns in order
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         cluster_distance_cols = sampleDists,
         col = colors,
         cluster_rows = FALSE,
         cluster_cols = FALSE)
#Heatmap with the rows and columns clustered
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         cluster_distance_cols = sampleDists,
         col = colors)

#Session info including packages used and their information
devtools::session_info()

#Heatmap of genes with most variance
topVarGenes <- head(order(rowVars(assay(vst_parametric)),decreasing=TRUE),20)
mat <- assay(vst_parametric)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(vst_parametric))
pheatmap(mat, annotation_col=df, labels_col = c("MD-E1", "MD-E2", "MD-E3", "MD-E4", "MD-E5", "MD-E6", "MD-E7", "MD-E8", 
                                                "Alc-F1", "Alc-F2", "Alc-F4", "Alc-F5", "Alc-F6", "Alc-F7", "Alc-F8", "Alc-F9"))
####Did Not Use Yet####
#Volcano Plot
# Make a basic volcano plot
#Genes that are highly dysregulated are farther to the left and right sides, 
#while highly significant changes appear higher on the plot.
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano Plot"))
# Add colored points: red if padj<0.05, orange of log2FC>1, green if both
with(subset(res, padj < 0.05), points(log2FoldChange, -log10(pvalue), pch=20, col="red")) 
#taking abs(absolute value) of log2 FC in order to see fold change greater than 1 in both negative and positive changes
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
# Label points with the textxy function from the calibrate plot
top <- head(sort(res$padj,decreasing=FALSE), n = 20)
max <- max(top)
with(subset(res, padj<max & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=resdata$GeneID, cex=.8))
