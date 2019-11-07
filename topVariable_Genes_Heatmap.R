# Get topvar genes of unique gene sets 
cluster1_25kb_rlogMatrix_Unique_topVarGenes_rownames <- head(cluster1_25kb_rlogMatrix_Unique_rownames[order(cluster1_25kb_rlogMatrix_Unique_rownames$RowVarAVG.x, decreasing = TRUE), ], nrow(cluster1_25kb_rlogMatrix_Unique_rownames)*.20)
cluster2_25kb_rlogMatrix_Unique_topVarGenes_rownames <- head(cluster2_25kb_rlogMatrix_Unique_rownames[order(cluster2_25kb_rlogMatrix_Unique_rownames$RowVarAVG.x, decreasing = TRUE), ], nrow(cluster2_25kb_rlogMatrix_Unique_rownames)*.20)
cluster3_25kb_rlogMatrix_Unique_topVarGenes_rownames <- head(cluster3_25kb_rlogMatrix_Unique_rownames[order(cluster3_25kb_rlogMatrix_Unique_rownames$RowVarAVG.x, decreasing = TRUE), ], nrow(cluster3_25kb_rlogMatrix_Unique_rownames)*.20)
cluster4_25kb_rlogMatrix_Unique_topVarGenes_rownames <- head(cluster4_25kb_rlogMatrix_Unique_rownames[order(cluster4_25kb_rlogMatrix_Unique_rownames$RowVarAVG.x, decreasing = TRUE), ], nrow(cluster4_25kb_rlogMatrix_Unique_rownames)*.20)
cluster5_25kb_rlogMatrix_Unique_topVarGenes_rownames <- head(cluster5_25kb_rlogMatrix_Unique_rownames[order(cluster5_25kb_rlogMatrix_Unique_rownames$RowVarAVG.x, decreasing = TRUE), ], nrow(cluster5_25kb_rlogMatrix_Unique_rownames)*.20)

# Center values around average for each row 
cluster1_25kb_rlogMatrix_Unique_topVarGenes_centered_rownames <- data.frame(cluster1_25kb_rlogMatrix_Unique_topVarGenes_rownames$Row.names, cluster1_25kb_rlogMatrix_Unique_topVarGenes_rownames[ ,2:11] - rowMeans(cluster1_25kb_rlogMatrix_Unique_topVarGenes_rownames[ ,2:11]))
cluster2_25kb_rlogMatrix_Unique_topVarGenes_centered_rownames <- data.frame(cluster2_25kb_rlogMatrix_Unique_topVarGenes_rownames$Row.names, cluster2_25kb_rlogMatrix_Unique_topVarGenes_rownames[ ,2:11] - rowMeans(cluster2_25kb_rlogMatrix_Unique_topVarGenes_rownames[ ,2:11]))
cluster3_25kb_rlogMatrix_Unique_topVarGenes_centered_rownames <- data.frame(cluster3_25kb_rlogMatrix_Unique_topVarGenes_rownames$Row.names, cluster3_25kb_rlogMatrix_Unique_topVarGenes_rownames[ ,2:11] - rowMeans(cluster3_25kb_rlogMatrix_Unique_topVarGenes_rownames[ ,2:11]))
cluster4_25kb_rlogMatrix_Unique_topVarGenes_centered_rownames <- data.frame(cluster4_25kb_rlogMatrix_Unique_topVarGenes_rownames$Row.names, cluster4_25kb_rlogMatrix_Unique_topVarGenes_rownames[ ,2:11] - rowMeans(cluster4_25kb_rlogMatrix_Unique_topVarGenes_rownames[ ,2:11]))
cluster5_25kb_rlogMatrix_Unique_topVarGenes_centered_rownames <- data.frame(cluster5_25kb_rlogMatrix_Unique_topVarGenes_rownames$Row.names, cluster5_25kb_rlogMatrix_Unique_topVarGenes_rownames[ ,2:11] - rowMeans(cluster5_25kb_rlogMatrix_Unique_topVarGenes_rownames[ ,2:11]))

# Fix column names 
colnames(cluster1_25kb_rlogMatrix_Unique_topVarGenes_centered_rownames) <- c("ENSGID", "Rep1_0", "Rep1_30", "Rep1_60", "Rep1_120", "Rep1_24", "Rep2_0", "Rep2_30", "Rep2_60", "Rep2_120", "Rep2_24")
colnames(cluster2_25kb_rlogMatrix_Unique_topVarGenes_centered_rownames) <- c("ENSGID", "Rep1_0", "Rep1_30", "Rep1_60", "Rep1_120", "Rep1_24", "Rep2_0", "Rep2_30", "Rep2_60", "Rep2_120", "Rep2_24")
colnames(cluster3_25kb_rlogMatrix_Unique_topVarGenes_centered_rownames) <- c("ENSGID", "Rep1_0", "Rep1_30", "Rep1_60", "Rep1_120", "Rep1_24", "Rep2_0", "Rep2_30", "Rep2_60", "Rep2_120", "Rep2_24")
colnames(cluster4_25kb_rlogMatrix_Unique_topVarGenes_centered_rownames) <- c("ENSGID", "Rep1_0", "Rep1_30", "Rep1_60", "Rep1_120", "Rep1_24", "Rep2_0", "Rep2_30", "Rep2_60", "Rep2_120", "Rep2_24")
colnames(cluster5_25kb_rlogMatrix_Unique_topVarGenes_centered_rownames) <- c("ENSGID", "Rep1_0", "Rep1_30", "Rep1_60", "Rep1_120", "Rep1_24", "Rep2_0", "Rep2_30", "Rep2_60", "Rep2_120", "Rep2_24")

# Get merged tables of peak annotations with meth information
Cluster1_peakAnnotation_methylation <- merge(cluster1_df, cluster1_methtable, by.y = "PeakID", by.x = "V4", all=FALSE)
Cluster2_peakAnnotation_methylation <- merge(cluster2_df, cluster2_methtable, by.y = "PeakID", by.x = "V4", all=FALSE)
Cluster3_peakAnnotation_methylation <- merge(cluster3_df, cluster3_methtable, by.y = "PeakID", by.x = "V4", all=FALSE)
Cluster4_peakAnnotation_methylation <- merge(cluster4_df, cluster4_methtable, by.y = "PeakID", by.x = "V4", all=FALSE)
Cluster5_peakAnnotation_methylation <- merge(cluster5_df, cluster5_methtable, by.y = "PeakID", by.x = "V4", all=FALSE)

# merge table above once more 
Cluster1_peakAnnotation_methylation_Unique_topVarGenes_Centered <- merge(Cluster1_peakAnnotation_methylation, cluster1_25kb_rlogMatrix_Unique_topVarGenes_centered_rownames, by.x = "ENSEMBL", by.y ="ENSGID", all=FALSE)
Cluster2_peakAnnotation_methylation_Unique_topVarGenes_Centered <- merge(Cluster2_peakAnnotation_methylation, cluster2_25kb_rlogMatrix_Unique_topVarGenes_centered_rownames, by.x = "ENSEMBL", by.y ="ENSGID", all=FALSE)
Cluster3_peakAnnotation_methylation_Unique_topVarGenes_Centered <- merge(Cluster3_peakAnnotation_methylation, cluster3_25kb_rlogMatrix_Unique_topVarGenes_centered_rownames, by.x = "ENSEMBL", by.y ="ENSGID", all=FALSE)
Cluster4_peakAnnotation_methylation_Unique_topVarGenes_Centered <- merge(Cluster4_peakAnnotation_methylation, cluster4_25kb_rlogMatrix_Unique_topVarGenes_centered_rownames, by.x = "ENSEMBL", by.y ="ENSGID", all=FALSE)
Cluster5_peakAnnotation_methylation_Unique_topVarGenes_Centered <- merge(Cluster5_peakAnnotation_methylation, cluster5_25kb_rlogMatrix_Unique_topVarGenes_centered_rownames, by.x = "ENSEMBL", by.y ="ENSGID", all=FALSE)

# Run pheatmap for clusters as a whole 

# also set up color palette breaks to center zero on white and maintain color scale for each cluster plot

paletteLength <- 50
myColor <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(paletteLength)
myBreaks1 <- c(seq(min(Cluster1_peakAnnotation_methylation_Unique_topVarGenes_Centered[,list]), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(Cluster1_peakAnnotation_methylation_Unique_topVarGenes_Centered[,list])/paletteLength, max(Cluster2_peakAnnotation_methylation_Unique_topVarGenes_Centered[,list]), length.out=floor(paletteLength/2)))
myBreaks2 <- c(seq(min(Cluster2_peakAnnotation_methylation_Unique_topVarGenes_Centered[,list]), 0, length.out=ceiling(paletteLength/2) + 1), 
               seq(max(Cluster2_peakAnnotation_methylation_Unique_topVarGenes_Centered[,list])/paletteLength, max(Cluster2_peakAnnotation_methylation_Unique_topVarGenes_Centered[,list]), length.out=floor(paletteLength/2)))
myBreaks3 <- c(seq(min(Cluster3_peakAnnotation_methylation_Unique_topVarGenes_Centered[,list]), 0, length.out=ceiling(paletteLength/2) + 1), 
               seq(max(Cluster3_peakAnnotation_methylation_Unique_topVarGenes_Centered[,list])/paletteLength, max(Cluster3_peakAnnotation_methylation_Unique_topVarGenes_Centered[,list]), length.out=floor(paletteLength/2)))
myBreaks4 <- c(seq(min(Cluster4_peakAnnotation_methylation_Unique_topVarGenes_Centered[,list]), 0, length.out=ceiling(paletteLength/2) + 1), 
               seq(max(Cluster4_peakAnnotation_methylation_Unique_topVarGenes_Centered[,list])/paletteLength, max(Cluster4_peakAnnotation_methylation_Unique_topVarGenes_Centered[,list]), length.out=floor(paletteLength/2)))
myBreaks5 <- c(seq(min(Cluster5_peakAnnotation_methylation_Unique_topVarGenes_Centered[,list]), 0, length.out=ceiling(paletteLength/2) + 1), 
               seq(max(Cluster5_peakAnnotation_methylation_Unique_topVarGenes_Centered[,list])/paletteLength, max(Cluster5_peakAnnotation_methylation_Unique_topVarGenes_Centered[,list]), length.out=floor(paletteLength/2)))

# Plot the heatmap
Cluster1_peakAnnotation_methylation_Unique_topVarGenes_Centered_All_Hclust <- pheatmap(Cluster1_peakAnnotation_methylation_Unique_topVarGenes_Centered[,list], cellwidth = 5, cellheight = 2.5, color = myColor, breaks = myBreaks1, cluster_cols=F, scale = "none", annotation_names_row=F, show_rownames=F, main = "ATAC Cluster 1 Nearby Genes <=25kb Away")
Cluster2_peakAnnotation_methylation_Unique_topVarGenes_Centered_All_Hclust <- pheatmap(Cluster2_peakAnnotation_methylation_Unique_topVarGenes_Centered[,list], cellwidth = 5, cellheight = 2.5, color = myColor, breaks = myBreaks2, cluster_cols=F, scale = "none", annotation_names_row=F, show_rownames=F, main = "ATAC Cluster 2 Nearby Genes <=25kb Away")
Cluster3_peakAnnotation_methylation_Unique_topVarGenes_Centered_All_Hclust <- pheatmap(Cluster3_peakAnnotation_methylation_Unique_topVarGenes_Centered[,list], cellwidth = 5, cellheight = 2.5, color = myColor, breaks = myBreaks3, cluster_cols=F, scale = "none", annotation_names_row=F, show_rownames=F, main = "ATAC Cluster 3 Nearby Genes <=25kb Away")
Cluster4_peakAnnotation_methylation_Unique_topVarGenes_Centered_All_Hclust <- pheatmap(Cluster4_peakAnnotation_methylation_Unique_topVarGenes_Centered[,list], cellwidth = 5, cellheight = 2.5, color = myColor, breaks = myBreaks4, cluster_cols=F, scale = "none", annotation_names_row=F, show_rownames=F, main = "ATAC Cluster 4 Nearby Genes <=25kb Away")
Cluster5_peakAnnotation_methylation_Unique_topVarGenes_Centered_All_Hclust <- pheatmap(Cluster5_peakAnnotation_methylation_Unique_topVarGenes_Centered[,list], cellwidth = 5, cellheight = 2.5, color = myColor, breaks = myBreaks5, cluster_cols=F, scale = "none", annotation_names_row=F, show_rownames=F, main = "ATAC Cluster 5 Nearby Genes <=25kb Away")

# Reorder original data the same as whole cluster heatmap
Cluster1_peakAnnotation_methylation_Unique_topVarGenes_Centered_ClusterOrdering <- Cluster1_peakAnnotation_methylation_Unique_topVarGenes_Centered[Cluster1_peakAnnotation_methylation_Unique_topVarGenes_Centered_All_Hclust$tree_row$order, ]
Cluster2_peakAnnotation_methylation_Unique_topVarGenes_Centered_ClusterOrdering <- Cluster2_peakAnnotation_methylation_Unique_topVarGenes_Centered[Cluster2_peakAnnotation_methylation_Unique_topVarGenes_Centered_All_Hclust$tree_row$order, ]
Cluster3_peakAnnotation_methylation_Unique_topVarGenes_Centered_ClusterOrdering <- Cluster3_peakAnnotation_methylation_Unique_topVarGenes_Centered[Cluster3_peakAnnotation_methylation_Unique_topVarGenes_Centered_All_Hclust$tree_row$order, ]
Cluster4_peakAnnotation_methylation_Unique_topVarGenes_Centered_ClusterOrdering <- Cluster4_peakAnnotation_methylation_Unique_topVarGenes_Centered[Cluster4_peakAnnotation_methylation_Unique_topVarGenes_Centered_All_Hclust$tree_row$order, ]
Cluster5_peakAnnotation_methylation_Unique_topVarGenes_Centered_ClusterOrdering <- Cluster5_peakAnnotation_methylation_Unique_topVarGenes_Centered[Cluster5_peakAnnotation_methylation_Unique_topVarGenes_Centered_All_Hclust$tree_row$order, ]

# Make heatmap of meth values without changing heatmap ordering from RNA Hclust
pheatmap(Cluster1_peakAnnotation_methylation_Unique_topVarGenes_Centered_ClusterOrdering$Meth, cellwidth = 5, cellheight = 2.5, color = colorRampPalette(rev(brewer.pal(n = 11, name = "PiYG")))(100), cluster_cols=F, cluster_rows = F, scale = "none", annotation_names_row=F, show_rownames=F, main = "ATAC Cluster 1 Methylation")
pheatmap(Cluster2_peakAnnotation_methylation_Unique_topVarGenes_Centered_ClusterOrdering$Meth, cellwidth = 5, cellheight = 2.5, color = colorRampPalette(rev(brewer.pal(n = 11, name = "PiYG")))(100), cluster_cols=F, cluster_rows = F, scale = "none", annotation_names_row=F, show_rownames=F, main = "ATAC Cluster 2 Methylation")
pheatmap(Cluster3_peakAnnotation_methylation_Unique_topVarGenes_Centered_ClusterOrdering$Meth, cellwidth = 5, cellheight = 2.5, color = colorRampPalette(rev(brewer.pal(n = 11, name = "PiYG")))(100), cluster_cols=F, cluster_rows = F, scale = "none", annotation_names_row=F, show_rownames=F, main = "ATAC Cluster 3 Methylation")
pheatmap(Cluster4_peakAnnotation_methylation_Unique_topVarGenes_Centered_ClusterOrdering$Meth, cellwidth = 5, cellheight = 2.5, color = colorRampPalette(rev(brewer.pal(n = 11, name = "PiYG")))(100), cluster_cols=F, cluster_rows = F, scale = "none", annotation_names_row=F, show_rownames=F, main = "ATAC Cluster 4 Methylation")
pheatmap(Cluster5_peakAnnotation_methylation_Unique_topVarGenes_Centered_ClusterOrdering$Meth, cellwidth = 5, cellheight = 2.5, color = colorRampPalette(rev(brewer.pal(n = 11, name = "PiYG")))(100), cluster_cols=F, cluster_rows = F, scale = "none", annotation_names_row=F, show_rownames=F, main = "ATAC Cluster 5 Methylation")
