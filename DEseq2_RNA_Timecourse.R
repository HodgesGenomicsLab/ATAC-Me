
# Load required R packages
library(DESeq2)
library(ashr)
library(dplyr)

# set directory location to where the count files from HTSseq are located
directory = '/Users/kellybarnett/Desktop/EH4889_noshift/DESeq2_RNA_Processing/count_files'
setwd(directory)

# Setup of experimental design and files
# Working directory should be directory of counts files from HTseq

# Get list of sample files
sampleFiles=grep('.IDcounts.txt',list.files(directory),value=TRUE)

# Get sample names
sampleNames=sub('.IDcounts.txt','',sampleFiles)

# Define time course conditions
sampleCondition <- factor( c( "0hr", "30min", "60min", "120min", "24hr",  "0hr", "30min", "60min", "120min", "24hr"))

# Build sample table with design info
sampleTable=data.frame(sampleName=sampleNames, fileName=sampleFiles, condition=sampleCondition)

# Build the DESeqDataSet from HTSeqCount table
cds = DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design= ~condition)

# Run standard differential expression analysis
# This is for most non-timecourse experimental design, see LRT for time course analysis
dds <- DESeq(cds)

# View results 
resultsNames(dds)

# Apply likelihood ratio test (LRT) to timecourse, this design assumes only time as a variant
ddsTC <- DESeq(dds, test="LRT", full = ~condition, reduced = ~1)

# Write results table of LRT test. Will contain Log2 fold change, p-values, and adjusted p-values for each gene
# Genes that have too little data to make a conclusion will be assigned a value of NA. This can be changed in the results() settings
resTC <- results(ddsTC)

# Log transforms the DESeq2 count table so that it is more appropriate for use in clustering or PCA
rld <- rlog(dds, blind=F) 

# Shrink the lfc of low count genes so they do not artificially skew results 
# These are shrunken LFC in reference to 0hr timepoint
resAsh30 <- lfcShrink(ddsTC, coef=4, type="ashr")
resAsh60 <- lfcShrink(ddsTC, coef=5, type="ashr")
resAsh120 <- lfcShrink(ddsTC, coef=2, type="ashr")
resAsh24 <- lfcShrink(ddsTC, coef=3, type="ashr")

# Build new dataframe of shrunken log fold changes
lfcShrink_df <- data.frame(rownames(ddsTC), as.numeric("0"), resAsh30$log2FoldChange, resAsh60$log2FoldChange, resAsh120$log2FoldChange, resAsh24$log2FoldChange, resTC$padj)
colnames(lfcShrink_df) <- c("ENSEMBL", "lfc0", "lfc30", "lfc60", "lfc120", "lfc24", "padj")

# Build similar dataframe but take the absolute value of log fold change
# Also add column that is the rowmax for absolute value lfc values 
lfcShrink_df_absvalue <- data.frame(rownames(ddsTC), as.numeric("0"), abs(resAsh30$log2FoldChange), abs(resAsh60$log2FoldChange), abs(resAsh120$log2FoldChange), abs(resAsh24$log2FoldChange), resTC$padj)
lfcShrink_df_absvalue <- data.frame(rownames(ddsTC), as.numeric("0"), abs(resAsh30$log2FoldChange), abs(resAsh60$log2FoldChange), abs(resAsh120$log2FoldChange), abs(resAsh24$log2FoldChange), apply(lfcShrink_df_absvalue[ ,2:6], 1, max), resTC$padj)
colnames(lfcShrink_df_absvalue) <- c("ENSEMBL", "lfc0", "lfc30", "lfc60", "lfc120", "lfc24", "lfcMax", "padj")

# Filter lfcShrink_df_absvalue for padj value to remove low confidence genes 
lfcShrink_df_absvalue <- filter(lfcShrink_df_absvalue, padj <= .005)

# Sort data by abs(logFC) maxima
lfcShrink_df_absvalue <- arrange(lfcShrink_df_absvalue, desc(lfcMax))

# Filter for logFC >=2 between any consecutive time points 
lfcShrink_df_absvalue_filt <- filter(lfcShrink_df_absvalue, lfcMax >= 2)


###### Filter data for LFC maxima 

# Filters the log transformed data against the lfc Maxima
rld.sigtest2 <- rld[ which(lfcShrink_df_absvalue$lfcMax >= 2), ]

# Generate matrix that we can feed to pheatmap for clustering
mat2 <- assay(rld.sigtest2)
# Remove rows that have no counts
mat2.1 <- mat2[ which(rowSums2(mat2) != 0), ]

# Re-order rows of  matrix for input into pheatmap
list <- c("KB1", "KB6", "KB2", "KB7", "KB3", "KB8", "KB4", "KB9", "KB5", "KBA")

# Kmeans clustering of log transformed counts data
# Optimal number of clusters for matrix of values was determined using gap statistic analysis
# Only genes that were beyond the p-adj filter and LFC filtering
heatmapres2.1<- pheatmap(mat2.1[,list], color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100), cellwidth = 30, cellheight = 30, cluster_cols=F, kmeans_k = 7, scale = "row", cutree_rows=7)


#  Generate heatmap with kmeans clustering and write image as SVG file
library(pheatmap)
library(RColorBrewer)
library(Cairo)
require(Cairo)
CairoSVG(file = "DESEQ_RNA_EH4889_7Clusters.svg", width = 7, height = 7)
datares <- pheatmap(mat2[,list], color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100), cellwidth = 30, cellheight = 30, cluster_cols=F, kmeans_k = 7, scale = "row", cutree_rows=7)
dev.off()

# Generate heatmap with hierarchical clustering and write image as SVG file
library(pheatmap)
library(RColorBrewer)
library(Cairo)
require(Cairo)
CairoSVG(file = "DESEQ_RNA_EH4889_HierCluster.svg", width = 7, height = 7)
datares_HC <- pheatmap(mat2[,list], color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100), cluster_cols=F, scale = "row", annotation_names_row=F)
dev.off()
