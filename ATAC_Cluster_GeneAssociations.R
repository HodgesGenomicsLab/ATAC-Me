# Load required libraries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(reshape2)


# Cluster 1 = 
# Cluster 2 = 
# Cluster 3 =
# Cluster 4 =
# Cluster 5 =

# Annotate bed files of genomic regions from TCseq k-meands clustering 
cluster1_PeakAnnotation <- annotatePeak("peak_subset_cluster_1.GREAT.bed", tssRegion = c(-2000, 500), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                         level = "transcript", assignGenomicAnnotation = TRUE,
                         genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                                                       "Downstream", "Intergenic"), annoDb = "org.Hs.eg.db", addFlankGeneInfo = TRUE,
                         flankDistance = 10000, sameStrand = FALSE, ignoreOverlap = FALSE,
                         ignoreUpstream = FALSE, ignoreDownstream = FALSE, overlap = "all",
                         verbose = TRUE)

cluster2_PeakAnnotation <- annotatePeak("peak_subset_cluster_2.GREAT.bed", tssRegion = c(-2000, 500), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                         level = "transcript", assignGenomicAnnotation = TRUE,
                         genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                                                       "Downstream", "Intergenic"), annoDb = "org.Hs.eg.db", addFlankGeneInfo = TRUE,
                         flankDistance = 10000, sameStrand = FALSE, ignoreOverlap = FALSE,
                         ignoreUpstream = FALSE, ignoreDownstream = FALSE, overlap = "all",
                         verbose = TRUE)

cluster3_PeakAnnotation <- annotatePeak("peak_subset_cluster_3.GREAT.bed", tssRegion = c(-2000, 500), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                         level = "transcript", assignGenomicAnnotation = TRUE,
                         genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                                                       "Downstream", "Intergenic"), annoDb = "org.Hs.eg.db", addFlankGeneInfo = TRUE,
                         flankDistance = 10000, sameStrand = FALSE, ignoreOverlap = FALSE,
                         ignoreUpstream = FALSE, ignoreDownstream = FALSE, overlap = "all",
                         verbose = TRUE)

cluster4_PeakAnnotation <- annotatePeak("peak_subset_cluster_4.GREAT.bed", tssRegion = c(-2000, 500), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                         level = "transcript", assignGenomicAnnotation = TRUE,
                         genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                                                       "Downstream", "Intergenic"), annoDb = "org.Hs.eg.db", addFlankGeneInfo = TRUE,
                         flankDistance = 10000, sameStrand = FALSE, ignoreOverlap = FALSE,
                         ignoreUpstream = FALSE, ignoreDownstream = FALSE, overlap = "all",
                         verbose = TRUE)

cluster5_PeakAnnotation <- annotatePeak("peak_subset_cluster_5.GREAT.bed", tssRegion = c(-2000, 500), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                         level = "transcript", assignGenomicAnnotation = TRUE,
                         genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                                                       "Downstream", "Intergenic"), annoDb = "org.Hs.eg.db", addFlankGeneInfo = TRUE,
                         flankDistance = 10000, sameStrand = FALSE, ignoreOverlap = FALSE,
                         ignoreUpstream = FALSE, ignoreDownstream = FALSE, overlap = "all",
                         verbose = TRUE)

all_dynamic_PeakAnnotation <- annotatePeak("all_cluster_peaks.bed", tssRegion = c(-2000, 500), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                            level = "transcript", assignGenomicAnnotation = TRUE,
                            genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                                                          "Downstream", "Intergenic"), annoDb = "org.Hs.eg.db", addFlankGeneInfo = TRUE,
                            flankDistance = 10000, sameStrand = FALSE, ignoreOverlap = FALSE,
                            ignoreUpstream = FALSE, ignoreDownstream = FALSE, overlap = "all",
                            verbose = TRUE)

all_static_PeakAnnotation <- annotatePeak("all_static_peaks.bed", tssRegion = c(-2000, 500), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                           level = "transcript", assignGenomicAnnotation = TRUE,
                           genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                                                         "Downstream", "Intergenic"), annoDb = "org.Hs.eg.db", addFlankGeneInfo = TRUE,
                           flankDistance = 10000, sameStrand = FALSE, ignoreOverlap = FALSE,
                           ignoreUpstream = FALSE, ignoreDownstream = FALSE, overlap = "all",
                           verbose = TRUE)

all_TCpeaks_PeakAnnotation <- annotatePeak("gf_tcseq.GREAT.bed", tssRegion = c(-2000, 500), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                            level = "transcript", assignGenomicAnnotation = TRUE,
                            genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                                                          "Downstream", "Intergenic"), annoDb = "org.Hs.eg.db", addFlankGeneInfo = TRUE,
                            flankDistance = 10000, sameStrand = FALSE, ignoreOverlap = FALSE,
                            ignoreUpstream = FALSE, ignoreDownstream = FALSE, overlap = "all",
                            verbose = TRUE)
# Genomic annotation barplots
cluster1_AnnoBar <- plotAnnoBar(cluster1_PeakAnnotation)
cluster2_AnnoBar <- plotAnnoBar(cluster2_PeakAnnotation)
cluster3_AnnoBar <- plotAnnoBar(cluster3_PeakAnnotation)
cluster4_AnnoBar <- plotAnnoBar(cluster4_PeakAnnotation)
cluster5_AnnoBar <- plotAnnoBar(cluster5_PeakAnnotation)
Static_AnnoBar <- plotAnnoBar(all_static_PeakAnnotation)
Dynamic_AnnoBar <- plotAnnoBar(all_dynamic_PeakAnnotation)
TCseqPeaks_AnnoBar <- plotAnnoBar(all_TCpeaks_PeakAnnotation)

# Convert the annotation table to an R data frame
# This data frame can be written as a text table 
cluster1_PeakAnnotation_df <- as.data.frame(cluster1_PeakAnnotation)
cluster2_PeakAnnotation_df <- as.data.frame(cluster2_PeakAnnotation)
cluster3_PeakAnnotation_df <- as.data.frame(cluster3_PeakAnnotation)
cluster4_PeakAnnotation_df <- as.data.frame(cluster4_PeakAnnotation)
cluster5_PeakAnnotation_df <- as.data.frame(cluster5_PeakAnnotation)
all_dynamic_PeakAnnotation_df <- as.data.frame(all_dynamic_PeakAnnotation)
all_static_PeakAnnotation_df <- as.data.frame(all_static_PeakAnnotation)
all_TCpeaks_PeakAnnotation_df <- as.data.frame(all_TCpeaks_PeakAnnotation)

# Filter for annotated nearby genes that are less than or equal to 25 kilobases away
cluster1_all_25kb_genes <- filter(cluster1_PeakAnnotation_df, abs(distanceToTSS) <= 25000)
cluster2_all_25kb_genes <- filter(cluster2_PeakAnnotation_df, abs(distanceToTSS) <= 25000)
cluster3_all_25kb_genes <- filter(cluster3_PeakAnnotation_df, abs(distanceToTSS) <= 25000)
cluster4_all_25kb_genes <- filter(cluster4_PeakAnnotation_df, abs(distanceToTSS) <= 25000)
cluster5_all_25kb_genes <- filter(cluster5_PeakAnnotation_df, abs(distanceToTSS) <= 25000)
all_dynamic_all_25kb_genes <- filter(all_dynamic_PeakAnnotation_df, abs(distanceToTSS) <= 25000)
all_static_all_25kb_genes <- filter(all_static_PeakAnnotation_df, abs(distanceToTSS) <= 25000)
all_TCpeaks_all_25kb_genes <- filter(all_TCpeaks_PeakAnnotation_df, abs(distanceToTSS) <= 25000)

# Add cluster group column to each table
cluster1_all_25kb_genes$Cluster <- c("cluster1")
cluster2_all_25kb_genes$Cluster <- c("cluster2")
cluster3_all_25kb_genes$Cluster <- c("cluster3")
cluster4_all_25kb_genes$Cluster <- c("cluster4")
cluster5_all_25kb_genes$Cluster <- c("cluster5")
all_dynamic_all_25kb_genes$Cluster <- c("all_dynamic")
all_static_all_25kb_genes$Cluster <- c("all_static")
all_TCpeaks_all_25kb_genes$Cluster <- c("all_TCpeaks")

# Combine table into one
ClusterAnnotation_25kb <- rbind(cluster1_all_25kb_genes, cluster2_all_25kb_genes, cluster3_all_25kb_genes, cluster4_all_25kb_genes, cluster5_all_25kb_genes, all_dynamic_all_25kb_genes, all_static_all_25kb_genes, all_TCpeaks_all_25kb_genes)

# Pull out counts for each ATAC peak cluster 
cluster1_25kb_rlogMatrix <- rld_matrix_df %>% filter(rownames(rld_matrix_df) %in% cluster1_all_25kb_genes$ENSEMBL)
cluster2_25kb_rlogMatrix <- rld_matrix_df %>% filter(rownames(rld_matrix_df) %in% cluster2_all_25kb_genes$ENSEMBL)
cluster3_25kb_rlogMatrix <- rld_matrix_df %>% filter(rownames(rld_matrix_df) %in% cluster3_all_25kb_genes$ENSEMBL)
cluster4_25kb_rlogMatrix <- rld_matrix_df %>% filter(rownames(rld_matrix_df) %in% cluster4_all_25kb_genes$ENSEMBL)
cluster5_25kb_rlogMatrix <- rld_matrix_df %>% filter(rownames(rld_matrix_df) %in% cluster5_all_25kb_genes$ENSEMBL)
all_dynamic_25kb_rlogMatrix <- rld_matrix_df %>% filter(rownames(rld_matrix_df) %in% all_dynamic_all_25kb_genes$ENSEMBL)
all_static_25kb_rlogMatrix <- rld_matrix_df %>% filter(rownames(rld_matrix_df) %in% all_static_all_25kb_genes$ENSEMBL)
all_TCpeaks_25kb_rlogMatrix <- rld_matrix_df %>% filter(rownames(rld_matrix_df) %in% all_TCpeaks_all_25kb_genes$ENSEMBL)

# Convert ENS ID back to rownames
cluster1_25kb_rlogMatrix <- data.frame(cluster1_25kb_rlogMatrix[,-1], row.names=cluster1_25kb_rlogMatrix[,1])
cluster2_25kb_rlogMatrix <- data.frame(cluster2_25kb_rlogMatrix[,-1], row.names=cluster2_25kb_rlogMatrix[,1])
cluster3_25kb_rlogMatrix <- data.frame(cluster3_25kb_rlogMatrix[,-1], row.names=cluster3_25kb_rlogMatrix[,1])
cluster4_25kb_rlogMatrix <- data.frame(cluster4_25kb_rlogMatrix[,-1], row.names=cluster4_25kb_rlogMatrix[,1])
cluster5_25kb_rlogMatrix <- data.frame(cluster5_25kb_rlogMatrix[,-1], row.names=cluster5_25kb_rlogMatrix[,1])
all_dynamic_25kb_rlogMatrix <- data.frame(all_dynamic_25kb_rlogMatrix[,-1], row.names=all_dynamic_25kb_rlogMatrix[,1])
all_static_25kb_rlogMatrix <- data.frame(all_static_25kb_rlogMatrix[,-1], row.names=all_static_25kb_rlogMatrix[,1])
all_TCpeaks_25kb_rlogMatrix <- data.frame(all_TCpeaks_25kb_rlogMatrix[,-1], row.names=all_TCpeaks_25kb_rlogMatrix[,1])

# Remove rows with all zeros from data.frame
cluster1_25kb_rlogMatrix <- cluster1_25kb_rlogMatrix[which(rowSums(cluster1_25kb_rlogMatrix) > 0),] 
cluster2_25kb_rlogMatrix <- cluster2_25kb_rlogMatrix[which(rowSums(cluster2_25kb_rlogMatrix) > 0),] 
cluster3_25kb_rlogMatrix <- cluster3_25kb_rlogMatrix[which(rowSums(cluster3_25kb_rlogMatrix) > 0),] 
cluster4_25kb_rlogMatrix <- cluster4_25kb_rlogMatrix[which(rowSums(cluster4_25kb_rlogMatrix) > 0),] 
cluster5_25kb_rlogMatrix <- cluster5_25kb_rlogMatrix[which(rowSums(cluster5_25kb_rlogMatrix) > 0),] 
all_dynamic_25kb_rlogMatrix <- all_dynamic_25kb_rlogMatrix[which(rowSums(all_dynamic_25kb_rlogMatrix) > 0),] 
all_static_25kb_rlogMatrix <- all_static_25kb_rlogMatrix[which(rowSums(all_static_25kb_rlogMatrix) > 0),] 
all_TCpeaks_25kb_rlogMatrix <- all_TCpeaks_25kb_rlogMatrix[which(rowSums(all_TCpeaks_25kb_rlogMatrix) > 0),] 

# Add row variance column to each table for replicates and also average variance of replicates
cluster1_25kb_rlogMatrix$RowVar1 <- rowVars(as.matrix(cluster1_25kb_rlogMatrix)[ , c("KB1", "KB2", "KB3", "KB4", "KB5")])
cluster2_25kb_rlogMatrix$RowVar1 <- rowVars(as.matrix(cluster2_25kb_rlogMatrix)[ , c("KB1", "KB2", "KB3", "KB4", "KB5")])
cluster3_25kb_rlogMatrix$RowVar1 <- rowVars(as.matrix(cluster3_25kb_rlogMatrix)[ , c("KB1", "KB2", "KB3", "KB4", "KB5")])
cluster4_25kb_rlogMatrix$RowVar1 <- rowVars(as.matrix(cluster4_25kb_rlogMatrix)[ , c("KB1", "KB2", "KB3", "KB4", "KB5")])
cluster5_25kb_rlogMatrix$RowVar1 <- rowVars(as.matrix(cluster5_25kb_rlogMatrix)[ , c("KB1", "KB2", "KB3", "KB4", "KB5")])
all_dynamic_25kb_rlogMatrix$RowVar1 <- rowVars(as.matrix(all_dynamic_25kb_rlogMatrix)[ , c("KB1", "KB2", "KB3", "KB4", "KB5")])
all_static_25kb_rlogMatrix$RowVar1 <- rowVars(as.matrix(all_static_25kb_rlogMatrix)[ , c("KB1", "KB2", "KB3", "KB4", "KB5")])
all_TCpeaks_25kb_rlogMatrix$RowVar1 <- rowVars(as.matrix(all_TCpeaks_25kb_rlogMatrix)[ , c("KB1", "KB2", "KB3", "KB4", "KB5")])

cluster1_25kb_rlogMatrix$RowVar2 <- rowVars(as.matrix(cluster1_25kb_rlogMatrix)[ , c("KB6", "KB7", "KB8", "KB9", "KBA")])
cluster2_25kb_rlogMatrix$RowVar2 <- rowVars(as.matrix(cluster2_25kb_rlogMatrix)[ , c("KB6", "KB7", "KB8", "KB9", "KBA")])
cluster3_25kb_rlogMatrix$RowVar2 <- rowVars(as.matrix(cluster3_25kb_rlogMatrix)[ , c("KB6", "KB7", "KB8", "KB9", "KBA")])
cluster4_25kb_rlogMatrix$RowVar2 <- rowVars(as.matrix(cluster4_25kb_rlogMatrix)[ , c("KB6", "KB7", "KB8", "KB9", "KBA")])
cluster5_25kb_rlogMatrix$RowVar2 <- rowVars(as.matrix(cluster5_25kb_rlogMatrix)[ , c("KB6", "KB7", "KB8", "KB9", "KBA")])
all_dynamic_25kb_rlogMatrix$RowVar2 <- rowVars(as.matrix(all_dynamic_25kb_rlogMatrix)[ , c("KB6", "KB7", "KB8", "KB9", "KBA")])
all_static_25kb_rlogMatrix$RowVar2 <- rowVars(as.matrix(all_static_25kb_rlogMatrix)[ , c("KB6", "KB7", "KB8", "KB9", "KBA")])
all_TCpeaks_25kb_rlogMatrix$RowVar2 <- rowVars(as.matrix(all_TCpeaks_25kb_rlogMatrix)[ , c("KB6", "KB7", "KB8", "KB9", "KBA")])

cluster1_25kb_rlogMatrix$RowVarAVG <- (cluster1_25kb_rlogMatrix$RowVar1+cluster1_25kb_rlogMatrix$RowVar2)/2
cluster2_25kb_rlogMatrix$RowVarAVG <- (cluster2_25kb_rlogMatrix$RowVar1+cluster2_25kb_rlogMatrix$RowVar2)/2
cluster3_25kb_rlogMatrix$RowVarAVG <- (cluster3_25kb_rlogMatrix$RowVar1+cluster3_25kb_rlogMatrix$RowVar2)/2
cluster4_25kb_rlogMatrix$RowVarAVG <- (cluster4_25kb_rlogMatrix$RowVar1+cluster4_25kb_rlogMatrix$RowVar2)/2
cluster5_25kb_rlogMatrix$RowVarAVG <- (cluster5_25kb_rlogMatrix$RowVar1+cluster5_25kb_rlogMatrix$RowVar2)/2
all_dynamic_25kb_rlogMatrix$RowVarAVG <- (all_dynamic_25kb_rlogMatrix$RowVar1+all_dynamic_25kb_rlogMatrix$RowVar2)/2
all_static_25kb_rlogMatrix$RowVarAVG <- (all_static_25kb_rlogMatrix$RowVar1+all_static_25kb_rlogMatrix$RowVar2)/2
all_TCpeaks_25kb_rlogMatrix$RowVarAVG <- (all_TCpeaks_25kb_rlogMatrix$RowVar1+all_TCpeaks_25kb_rlogMatrix$RowVar2)/2

# order lists according to row variance and take top 20% most variable genes
cluster1_25kb_rlogMatrix_topVarGenes <- head(cluster1_25kb_rlogMatrix[order(cluster1_25kb_rlogMatrix$RowVarAVG, decreasing = TRUE), ], nrow(cluster1_25kb_rlogMatrix)*.20)
cluster2_25kb_rlogMatrix_topVarGenes <- head(cluster2_25kb_rlogMatrix[order(cluster2_25kb_rlogMatrix$RowVarAVG, decreasing = TRUE), ], nrow(cluster2_25kb_rlogMatrix)*.20)
cluster3_25kb_rlogMatrix_topVarGenes <- head(cluster3_25kb_rlogMatrix[order(cluster3_25kb_rlogMatrix$RowVarAVG, decreasing = TRUE), ], nrow(cluster3_25kb_rlogMatrix)*.20)
cluster4_25kb_rlogMatrix_topVarGenes <- head(cluster4_25kb_rlogMatrix[order(cluster4_25kb_rlogMatrix$RowVarAVG, decreasing = TRUE), ], nrow(cluster4_25kb_rlogMatrix)*.20)
cluster5_25kb_rlogMatrix_topVarGenes <- head(cluster5_25kb_rlogMatrix[order(cluster5_25kb_rlogMatrix$RowVarAVG, decreasing = TRUE), ], nrow(cluster5_25kb_rlogMatrix)*.20)
all_dynamic_25kb_rlogMatrix_topVarGenes <- head(all_dynamic_25kb_rlogMatrix[order(all_dynamic_25kb_rlogMatrix$RowVarAVG, decreasing = TRUE), ], nrow(all_dynamic_25kb_rlogMatrix)*.20)
all_static_25kb_rlogMatrix_topVarGenes <- head(all_static_25kb_rlogMatrix[order(all_static_25kb_rlogMatrix$RowVarAVG, decreasing = TRUE), ], nrow(all_static_25kb_rlogMatrix)*.20)
all_TCpeaks_25kb_rlogMatrix_topVarGenes <- head(all_TCpeaks_25kb_rlogMatrix[order(all_TCpeaks_25kb_rlogMatrix$RowVarAVG, decreasing = TRUE), ], nrow(all_TCpeaks_25kb_rlogMatrix)*.20)

# Center matrix values around average of each row
cluster1_25kb_rlogMatrix_topVarGenes_centered <- cluster1_25kb_rlogMatrix_topVarGenes[ ,1:10] - rowMeans(cluster1_25kb_rlogMatrix_topVarGenes[ ,1:10])
cluster2_25kb_rlogMatrix_topVarGenes_centered <- cluster2_25kb_rlogMatrix_topVarGenes[ ,1:10] - rowMeans(cluster2_25kb_rlogMatrix_topVarGenes[ ,1:10])
cluster3_25kb_rlogMatrix_topVarGenes_centered <- cluster3_25kb_rlogMatrix_topVarGenes[ ,1:10] - rowMeans(cluster3_25kb_rlogMatrix_topVarGenes[ ,1:10])
cluster4_25kb_rlogMatrix_topVarGenes_centered <- cluster4_25kb_rlogMatrix_topVarGenes[ ,1:10] - rowMeans(cluster4_25kb_rlogMatrix_topVarGenes[ ,1:10])
cluster5_25kb_rlogMatrix_topVarGenes_centered <- cluster5_25kb_rlogMatrix_topVarGenes[ ,1:10] - rowMeans(cluster5_25kb_rlogMatrix_topVarGenes[ ,1:10])
all_dynamic_25kb_rlogMatrix_topVarGenes_centered <- all_dynamic_25kb_rlogMatrix_topVarGenes[ ,1:10] - rowMeans(all_dynamic_25kb_rlogMatrix_topVarGenes[ ,1:10])
all_static_25kb_rlogMatrix_topVarGenes_centered <- all_static_25kb_rlogMatrix_topVarGenes[ ,1:10] - rowMeans(all_static_25kb_rlogMatrix_topVarGenes[ ,1:10])
all_TCpeaks_25kb_rlogMatrix_topVarGenes_centered <- all_TCpeaks_25kb_rlogMatrix_topVarGenes[ ,1:10] - rowMeans(all_TCpeaks_25kb_rlogMatrix_topVarGenes[ ,1:10])

#### Get heatmaps of unique gene lists 

# Get set of genes that occur in only 1 cluster
All_Cluster_Genes_withVariance <- rbind(cluster1_25kb_rlogMatrix, cluster2_25kb_rlogMatrix, cluster3_25kb_rlogMatrix, cluster4_25kb_rlogMatrix, cluster5_25kb_rlogMatrix)
All_Cluster_Genes_withVariance_UniqueOnly <- All_Cluster_Genes_withVariance[!(duplicated(All_Cluster_Genes_withVariance) | duplicated(All_Cluster_Genes_withVariance, fromLast = TRUE)), ]

# Get rid of non-unique genes
cluster1_25kb_rlogMatrix_Unique <- intersect(cluster1_25kb_rlogMatrix, All_Cluster_Genes_withVariance_UniqueOnly)
cluster2_25kb_rlogMatrix_Unique <- intersect(cluster2_25kb_rlogMatrix, All_Cluster_Genes_withVariance_UniqueOnly)
cluster3_25kb_rlogMatrix_Unique <- intersect(cluster3_25kb_rlogMatrix, All_Cluster_Genes_withVariance_UniqueOnly)
cluster4_25kb_rlogMatrix_Unique <- intersect(cluster4_25kb_rlogMatrix, All_Cluster_Genes_withVariance_UniqueOnly)
cluster5_25kb_rlogMatrix_Unique <- intersect(cluster5_25kb_rlogMatrix, All_Cluster_Genes_withVariance_UniqueOnly)
all_dynamic_25kb_rlogMatrix_Unique <- intersect(all_dynamic_25kb_rlogMatrix, All_Cluster_Genes_withVariance_UniqueOnly)
all_static_25kb_rlogMatrix_Unique <- intersect(all_static_25kb_rlogMatrix, All_Cluster_Genes_withVariance_UniqueOnly)
all_TCpeaks_25kb_rlogMatrix_Unique <- intersect(all_TCpeaks_25kb_rlogMatrix, All_Cluster_Genes_withVariance_UniqueOnly)

# Preserve rownames to make downstream analysis easier
cluster1_25kb_rlogMatrix_Unique_rownames <- merge(cluster1_25kb_rlogMatrix, All_Cluster_Genes_withVariance_UniqueOnly, by.x = 0, by.y = 0, all=FALSE)[ ,1:14]
cluster2_25kb_rlogMatrix_Unique_rownames <- merge(cluster2_25kb_rlogMatrix, All_Cluster_Genes_withVariance_UniqueOnly, by.x = 0, by.y = 0, all=FALSE)[ ,1:14]
cluster3_25kb_rlogMatrix_Unique_rownames <- merge(cluster3_25kb_rlogMatrix, All_Cluster_Genes_withVariance_UniqueOnly, by.x = 0, by.y = 0, all=FALSE)[ ,1:14]
cluster4_25kb_rlogMatrix_Unique_rownames <- merge(cluster4_25kb_rlogMatrix, All_Cluster_Genes_withVariance_UniqueOnly, by.x = 0, by.y = 0, all=FALSE)[ ,1:14]
cluster5_25kb_rlogMatrix_Unique_rownames <- merge(cluster5_25kb_rlogMatrix, All_Cluster_Genes_withVariance_UniqueOnly, by.x = 0, by.y = 0, all=FALSE)[ ,1:14]
all_dynamic_25kb_rlogMatrix_Unique_rownames <- merge(all_dynamic_25kb_rlogMatrix, All_Cluster_Genes_withVariance_UniqueOnly, by.x = 0, by.y = 0, all=FALSE)[ ,1:14]
all_static_25kb_rlogMatrix_Unique_rownames <- merge(all_static_25kb_rlogMatrix, All_Cluster_Genes_withVariance_UniqueOnly, by.x = 0, by.y = 0, all=FALSE)[ ,1:14]
all_TCpeaks_25kb_rlogMatrix_Unique_rownames <- merge(all_TCpeaks_25kb_rlogMatrix, All_Cluster_Genes_withVariance_UniqueOnly, by.x = 0, by.y = 0, all=FALSE)[ ,1:14]


# Get topvar genes of unique gene sets 

cluster1_25kb_rlogMatrix_Unique_topVarGenes <- head(cluster1_25kb_rlogMatrix_Unique[order(cluster1_25kb_rlogMatrix_Unique$RowVarAVG, decreasing = TRUE), ], nrow(cluster1_25kb_rlogMatrix_Unique)*.20)
cluster2_25kb_rlogMatrix_Unique_topVarGenes <- head(cluster2_25kb_rlogMatrix_Unique[order(cluster2_25kb_rlogMatrix_Unique$RowVarAVG, decreasing = TRUE), ], nrow(cluster2_25kb_rlogMatrix_Unique)*.20)
cluster3_25kb_rlogMatrix_Unique_topVarGenes <- head(cluster3_25kb_rlogMatrix_Unique[order(cluster3_25kb_rlogMatrix_Unique$RowVarAVG, decreasing = TRUE), ], nrow(cluster3_25kb_rlogMatrix_Unique)*.20)
cluster4_25kb_rlogMatrix_Unique_topVarGenes <- head(cluster4_25kb_rlogMatrix_Unique[order(cluster4_25kb_rlogMatrix_Unique$RowVarAVG, decreasing = TRUE), ], nrow(cluster4_25kb_rlogMatrix_Unique)*.20)
cluster5_25kb_rlogMatrix_Unique_topVarGenes <- head(cluster5_25kb_rlogMatrix_Unique[order(cluster5_25kb_rlogMatrix_Unique$RowVarAVG, decreasing = TRUE), ], nrow(cluster5_25kb_rlogMatrix_Unique)*.20)
all_dynamic_25kb_rlogMatrix_Unique_topVarGenes <- head(all_dynamic_25kb_rlogMatrix_Unique[order(all_dynamic_25kb_rlogMatrix_Unique$RowVarAVG, decreasing = TRUE), ], nrow(all_dynamic_25kb_rlogMatrix_Unique)*.20)
all_static_25kb_rlogMatrix_Unique_topVarGenes <- head(all_static_25kb_rlogMatrix_Unique[order(all_static_25kb_rlogMatrix_Unique$RowVarAVG, decreasing = TRUE), ], nrow(all_static_25kb_rlogMatrix_Unique)*.20)
all_TCpeaks_25kb_rlogMatrix_Unique_topVarGenes <- head(all_TCpeaks_25kb_rlogMatrix_Unique[order(all_TCpeaks_25kb_rlogMatrix_Unique$RowVarAVG, decreasing = TRUE), ], nrow(all_TCpeaks_25kb_rlogMatrix_Unique)*.20)

# Center values around average for each row 
cluster1_25kb_rlogMatrix_Unique_topVarGenes_centered <- cluster1_25kb_rlogMatrix_Unique_topVarGenes[ ,1:10] - rowMeans(cluster1_25kb_rlogMatrix_Unique_topVarGenes[ ,1:10])
cluster2_25kb_rlogMatrix_Unique_topVarGenes_centered <- cluster2_25kb_rlogMatrix_Unique_topVarGenes[ ,1:10] - rowMeans(cluster2_25kb_rlogMatrix_Unique_topVarGenes[ ,1:10])
cluster3_25kb_rlogMatrix_Unique_topVarGenes_centered <- cluster3_25kb_rlogMatrix_Unique_topVarGenes[ ,1:10] - rowMeans(cluster3_25kb_rlogMatrix_Unique_topVarGenes[ ,1:10])
cluster4_25kb_rlogMatrix_Unique_topVarGenes_centered <- cluster4_25kb_rlogMatrix_Unique_topVarGenes[ ,1:10] - rowMeans(cluster4_25kb_rlogMatrix_Unique_topVarGenes[ ,1:10])
cluster5_25kb_rlogMatrix_Unique_topVarGenes_centered <- cluster5_25kb_rlogMatrix_Unique_topVarGenes[ ,1:10] - rowMeans(cluster5_25kb_rlogMatrix_Unique_topVarGenes[ ,1:10])
all_dynamic_25kb_rlogMatrix_Unique_topVarGenes_centered <- all_dynamic_25kb_rlogMatrix_Unique_topVarGenes[ ,1:10] - rowMeans(all_dynamic_25kb_rlogMatrix_Unique_topVarGenes[ ,1:10])
all_static_25kb_rlogMatrix_Unique_topVarGenes_centered <- all_static_25kb_rlogMatrix_Unique_topVarGenes[ ,1:10] - rowMeans(all_static_25kb_rlogMatrix_Unique_topVarGenes[ ,1:10])
all_TCpeaks_25kb_rlogMatrix_Unique_topVarGenes_centered <- all_TCpeaks_25kb_rlogMatrix_Unique_topVarGenes[ ,1:10] - rowMeans(all_TCpeaks_25kb_rlogMatrix_Unique_topVarGenes[ ,1:10])

# Generate pheatmap of topvar genes of unique cluster gene sets
# Hierarchical clustering 
cluster1_25kb_rlogMatrix_Unique_HClust_topVar <- pheatmap(cluster1_25kb_rlogMatrix_Unique_topVarGenes_centered[,list], color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100), cluster_cols=F, scale = "none", annotation_names_row=F, show_rownames=F, main = "ATAC Cluster 1 Nearby Genes <=25kb Away")
cluster2_25kb_rlogMatrix_Unique_HClust_topVar <- pheatmap(cluster2_25kb_rlogMatrix_Unique_topVarGenes_centered[,list], color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100), cluster_cols=F, scale = "none", annotation_names_row=F, show_rownames=F, main = "ATAC Cluster 2 Nearby Genes <=25kb Away")
cluster3_25kb_rlogMatrix_Unique_HClust_topVar <- pheatmap(cluster3_25kb_rlogMatrix_Unique_topVarGenes_centered[,list], color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100), cluster_cols=F, scale = "none", annotation_names_row=F, show_rownames=F, main = "ATAC Cluster 3 Nearby Genes <=25kb Away")
cluster4_25kb_rlogMatrix_Unique_HClust_topVar <- pheatmap(cluster4_25kb_rlogMatrix_Unique_topVarGenes_centered[,list], color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100), cluster_cols=F, scale = "none", annotation_names_row=F, show_rownames=F, main = "ATAC Cluster 4 Nearby Genes <=25kb Away")
cluster5_25kb_rlogMatrix_Unique_HClust_topVar <- pheatmap(cluster5_25kb_rlogMatrix_Unique_topVarGenes_centered[,list], color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100), cluster_cols=F, scale = "none", annotation_names_row=F, show_rownames=F, main = "ATAC Cluster 5 Nearby Genes <=25kb Away")
all_dynamic_25kb_rlogMatrix_Unique_HClust_topVar <- pheatmap(all_dynamic_25kb_rlogMatrix_Unique_topVarGenes_centered[,list], color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100), cluster_cols=F, scale = "none", annotation_names_row=F, show_rownames=F, main = "ATAC All Dynamic Nearby Genes <=25kb Away")
all_static_25kb_rlogMatrix_Unique_HClust_topVar <- pheatmap(all_static_25kb_rlogMatrix_Unique_topVarGenes_centered[,list], color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100), cluster_cols=F, scale = "none", annotation_names_row=F, show_rownames=F, main = "ATAC All Static Nearby Genes <=25kb Away")
all_TCpeaks_25kb_rlogMatrix_Unique_HClust_topVar <- pheatmap(all_TCpeaks_25kb_rlogMatrix_Unique_topVarGenes_centered[,list], color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100), cluster_cols=F, scale = "none", annotation_names_row=F, show_rownames=F, main = "ATAC All TCseq Nearby Genes <=25kb Away")
