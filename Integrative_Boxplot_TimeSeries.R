library(dplyr)
library(ggplot2)

######### Read in and prep files ######### 

# Directory of meth files 
# /data/hodges_lab/EH4889/fastqs/trimmed_reads/temp/mapped_reads/methprocessed
# Relevant Files in above directory:
ATAC-Me_THP1_PMA_MERGE_0.meth
ATAC-Me_THP1_PMA_MERGE_30.meth
ATAC-Me_THP1_PMA_MERGE_60.meth
ATAC-Me_THP1_PMA_MERGE_120.meth
ATAC-Me_THP1_PMA_MERGE_24hr.meth
# Make above files into bed files so that we can intersect them with cluster regions using bedtools
awk '{OFS="\t"; print $1,$2,$2+1,$5,$6}' ATAC-Me_THP1_PMA_MERGE_0.meth > ATAC-Me_THP1_PMA_MERGE_0.meth.bed 
awk '{OFS="\t"; print $1,$2,$2+1,$5,$6}' ATAC-Me_THP1_PMA_MERGE_30.meth > ATAC-Me_THP1_PMA_MERGE_30.meth.bed 
awk '{OFS="\t"; print $1,$2,$2+1,$5,$6}' ATAC-Me_THP1_PMA_MERGE_60.meth > ATAC-Me_THP1_PMA_MERGE_60.meth.bed 
awk '{OFS="\t"; print $1,$2,$2+1,$5,$6}' ATAC-Me_THP1_PMA_MERGE_120.meth > ATAC-Me_THP1_PMA_MERGE_120.meth.bed 
awk '{OFS="\t"; print $1,$2,$2+1,$5,$6}' ATAC-Me_THP1_PMA_MERGE_24hr.meth > ATAC-Me_THP1_PMA_MERGE_24hr.meth.bed 

# Directory of cluster bed files 
# /data/hodges_lab/EH4889/trimmed_reads/mapped_reads/methprocessed/methylKit_files/TC_seq_EH/
# Relevant files in above directory:
peak_subset_cluster_1.GREAT.bed  
peak_subset_cluster_2.GREAT.bed  
peak_subset_cluster_3.GREAT.bed  
peak_subset_cluster_4.GREAT.bed  
peak_subset_cluster_5.GREAT.bed 
all_cluster_peaks.bed
all_static_peaks.bed
gf_tcseq.GREAT.bed

# Get CpGs that are in clusters using bedtools intersect
bedtools intersect -wo -a ATAC-Me_THP1_PMA_MERGE_0.meth.bed -b /data/hodges_lab/EH4889/trimmed_reads/mapped_reads/methprocessed/methylKit_files/TC_seq_EH/all_cluster_peaks.bed > all_cluster_peaks_CpGs_0_meth.bed
bedtools intersect -wo -a ATAC-Me_THP1_PMA_MERGE_30.meth.bed -b /data/hodges_lab/EH4889/trimmed_reads/mapped_reads/methprocessed/methylKit_files/TC_seq_EH/all_cluster_peaks.bed > all_cluster_peaks_CpGs_30_meth.bed
bedtools intersect -wo -a ATAC-Me_THP1_PMA_MERGE_60.meth.bed -b /data/hodges_lab/EH4889/trimmed_reads/mapped_reads/methprocessed/methylKit_files/TC_seq_EH/all_cluster_peaks.bed > all_cluster_peaks_CpGs_60_meth.bed
bedtools intersect -wo -a ATAC-Me_THP1_PMA_MERGE_120.meth.bed -b /data/hodges_lab/EH4889/trimmed_reads/mapped_reads/methprocessed/methylKit_files/TC_seq_EH/all_cluster_peaks.bed > all_cluster_peaks_CpGs_120_meth.bed
bedtools intersect -wo -a ATAC-Me_THP1_PMA_MERGE_24hr.meth.bed -b /data/hodges_lab/EH4889/trimmed_reads/mapped_reads/methprocessed/methylKit_files/TC_seq_EH/all_cluster_peaks.bed > all_cluster_peaks_CpGs_24hr_meth.bed

# Read tables into R 
setwd("/Users/kellybarnett/Desktop/EH4889_noshift/")
all_cluster_peaks_CpGs_0_meth <- read.delim2(file = "all_cluster_peaks_CpGs_0_meth.bed", header = FALSE)
all_cluster_peaks_CpGs_30_meth <- read.delim2(file = "all_cluster_peaks_CpGs_30_meth.bed", header = FALSE)
all_cluster_peaks_CpGs_60_meth <- read.delim2(file = "all_cluster_peaks_CpGs_60_meth.bed", header = FALSE)
all_cluster_peaks_CpGs_120_meth <- read.delim2(file = "all_cluster_peaks_CpGs_120_meth.bed", header = FALSE)
all_cluster_peaks_CpGs_24hr_meth <- read.delim2(file = "all_cluster_peaks_CpGs_24hr_meth.bed", header = FALSE)

# Add time column
all_cluster_peaks_CpGs_0_meth$Time <- 0
all_cluster_peaks_CpGs_30_meth$Time <- 0.5
all_cluster_peaks_CpGs_60_meth$Time <- 1
all_cluster_peaks_CpGs_120_meth$Time <- 2
all_cluster_peaks_CpGs_24hr_meth$Time <- 24


######### Make accessibility table and standardize data ######### 

# Center values around mean and standardize
TCseq_Clustering_Counts_Table <- data.frame(TCseq_Clustering_Counts_Table, apply(TCseq_Clustering_Counts_Table[,2:6], 1, sd))
colnames(TCseq_Clustering_Counts_Table) <- c("PeakID", "0", "0.5", "1", "2", "24", "Chr", "Start", "End", "Cluster", "StDev")
TCseq_Clustering_Counts_Table_Standardized <- data.frame(TCseq_Clustering_Counts_Table[,1], (TCseq_Clustering_Counts_Table[,2:6] - rowMeans(TCseq_Clustering_Counts_Table[,2:6]))/(TCseq_Clustering_Counts_Table$StDev), TCseq_Clustering_Counts_Table[,7:10])
# TCseq_Clustering_Counts_Table_Centered <- data.frame(TCseq_Clustering_Counts_Table[,1], TCseq_Clustering_Counts_Table[,2:6] - rowMeans(TCseq_Clustering_Counts_Table[,2:6])/apply(TCseq_Clustering_Counts_Table[,2:6], 1, sd), TCseq_Clustering_Counts_Table[,7:10])
# Add column names
colnames(TCseq_Clustering_Counts_Table_Standardized) <- c("PeakID", "0", "0.5", "1", "2", "24", "Chr", "Start", "End", "Cluster")
# Change to tidy format table
#gather(TCseq_Clustering_Counts_Table_Centered, "Time", "ATAC_reads", 2:6) -> TCseq_Clustering_Counts_Table_tidy
gather(TCseq_Clustering_Counts_Table_Standardized, "Time", "ATAC_reads", 2:6) -> TCseq_Clustering_Counts_Table_Standardized_tidy

######### Scatterplot of Meth versus Accessibility ######### 
# Merge ATAC_counts table with meth table
ATAC_and_Meth_ClusterCpG_Table <- merge(meth_time_table, TCseq_Clustering_Counts_Table_tidy, by = c("PeakID", "Chr", "Start", "End", "Cluster", "Time"))
# Change meth column from factor type to numeric
ATAC_and_Meth_ClusterCpG_Table$Meth <- as.numeric(as.character(ATAC_and_Meth_ClusterCpG_Table$Meth))

### Generate plot for Cluster CpGs
Clust1_Plot <- ggplot(filter(ATAC_and_Meth_ClusterCpG_Table, Cluster == 3 & CpG_Reads >= 5), aes(x=Meth, y=ATAC_reads))+
  ggtitle("Meth vs. Access")+
  geom_point(alpha=0.025, color="steelblue")+
  ylab("ATAC Read Count Deviation from Mean")+
  xlab("CpG Methylation")+
  ylim(-50, 75)+
  xlim(0, 1)+
  theme(plot.title = element_text(lineheight=.8, face="bold", vjust=1),
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust=0.3))+ theme_bw()

# Facet wrap according to time for cluster
Clust1_Plot + facet_wrap( ~ Time, scales = "fixed")

########## Standardize methylation data ########## 

# Bind tables together to make tidy
meth_time_table <- rbind(all_cluster_peaks_CpGs_0_meth, all_cluster_peaks_CpGs_30_meth, all_cluster_peaks_CpGs_60_meth, all_cluster_peaks_CpGs_120_meth, all_cluster_peaks_CpGs_24hr_meth)
# Add column names
colnames(meth_time_table) <- c("CpG_Chr", "CpG_Start", "CpG_End", "Meth", "CpG_Reads", "Chr", "Start", "End", "PeakID","Cluster", "overlap", "Time")

# Spread meth_time_table to wide format that we can easily center the data in each row
meth_time_table_wideformat <- data.frame(all_cluster_peaks_CpGs_0_meth[,c(1,2,3,5,7,8.9,10)] ,all_cluster_peaks_CpGs_0_meth$V4, all_cluster_peaks_CpGs_30_meth$V4, all_cluster_peaks_CpGs_60_meth$V4, all_cluster_peaks_CpGs_120_meth$V4, all_cluster_peaks_CpGs_24hr_meth$V4)
# add column names
colnames(meth_time_table_wideformat) <- c("CpG_Chr", "CpG_Start", "CpG_End", "CpG_Reads", "Start", "End", "Cluster", "0", "0.5", "1", "2", "24")
# Change time column to numeric data type 
meth_time_table_wideformat$`0` <- as.numeric(as.character(meth_time_table_wideformat$`0`))
meth_time_table_wideformat$`0.5` <- as.numeric(as.character(meth_time_table_wideformat$`0.5`))
meth_time_table_wideformat$`1` <- as.numeric(as.character(meth_time_table_wideformat$`1`))
meth_time_table_wideformat$`2` <- as.numeric(as.character(meth_time_table_wideformat$`2`))
meth_time_table_wideformat$`24` <- as.numeric(as.character(meth_time_table_wideformat$`24`))
# Center and standardize the methylation levels 
meth_time_table_wideformat <- data.frame(meth_time_table_wideformat, apply(meth_time_table_wideformat[,8:12], 1, sd))
colnames(meth_time_table_wideformat) <- c("CpG_Chr", "CpG_Start", "CpG_End", "CpG_Reads", "Start", "End", "Cluster", "0", "0.5", "1", "2", "24", "StDev")
# meth_time_table_wideformat_standardized <- (meth_time_table_wideformat[,8:12] - meth_time_table_wideformat[,8])/(meth_time_table_wideformat$StDev)
meth_time_table_wideformat_standardized <- (meth_time_table_wideformat[,8:12] - rowMeans(meth_time_table_wideformat[,8:12]))/(meth_time_table_wideformat$StDev)
meth_time_table_wideformat_standardized <- data.frame(meth_time_table_wideformat[,1:7], meth_time_table_wideformat_standardized)
colnames(meth_time_table_wideformat_standardized) <- c("CpG_Chr", "CpG_Start", "CpG_End", "CpG_Reads", "Start", "End", "Cluster", "0", "0.5", "1", "2", "24") 

# Define function to filter our nan values due to divide by zero
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
# Apply filter to table
meth_time_table_wideformat_standardized[is.nan(meth_time_table_wideformat_standardized)] <- 0
# Filter out low count CpG_Reads
meth_time_table_wideformat_standardized_filter5 <- filter(meth_time_table_wideformat_standardized, CpG_Reads >= 5)
# Turn table into long format for easier  plotting
meth_time_table_centered_longformat <- gather(meth_time_table_wideformat_standardized, "Time", "Meth_stand", 8:12)
meth_time_table_centered_longformat_filter5 <- gather(meth_time_table_wideformat_standardized_filter5, "Time", "Meth_stand", 8:12)

# Change Time to numeric data type
meth_time_table_centered_longformat$Time <- as.factor(meth_time_table_centered_longformat$Time)
meth_time_table_centered_longformat_filter5$Time <- as.factor(meth_time_table_centered_longformat_filter5$Time)

######### make standardized table of nearby gene expression levels #########  

# Add cluster column
cluster1_25kb_rlogMatrix_Unique_topVarGenes_rownames$Cluster <- 1
cluster2_25kb_rlogMatrix_Unique_topVarGenes_rownames$Cluster <- 2
cluster3_25kb_rlogMatrix_Unique_topVarGenes_rownames$Cluster <- 3
cluster4_25kb_rlogMatrix_Unique_topVarGenes_rownames$Cluster <- 4
cluster5_25kb_rlogMatrix_Unique_topVarGenes_rownames$Cluster <- 5

# Bind all clusters together in a single table
All_clusters_25kb_rlogMatrix_Unique_topVarGenes_rownames <- rbind(cluster1_25kb_rlogMatrix_Unique_topVarGenes_rownames,
      cluster2_25kb_rlogMatrix_Unique_topVarGenes_rownames,
      cluster3_25kb_rlogMatrix_Unique_topVarGenes_rownames,
      cluster4_25kb_rlogMatrix_Unique_topVarGenes_rownames,
      cluster5_25kb_rlogMatrix_Unique_topVarGenes_rownames)

# Fix column names 
colnames(All_clusters_25kb_rlogMatrix_Unique_topVarGenes_rownames) <- c("ENSGID", "Rep1_0", "Rep1_30", "Rep1_60", "Rep1_120", "Rep1_24", "Rep2_0", "Rep2_30", "Rep2_60", "Rep2_120", "Rep2_24", "RowVar1", "RowVar2", "RowVarAVG", "Cluster")

# Make a table of replicate averages 

All_clusters_25kb_rlogMatrix_Unique_topVarGenes_Avg <- data.frame(All_clusters_25kb_rlogMatrix_Unique_topVarGenes_rownames[,1],
           (All_clusters_25kb_rlogMatrix_Unique_topVarGenes_rownames[,2] + All_clusters_25kb_rlogMatrix_Unique_topVarGenes_rownames[,7])/2,
           (All_clusters_25kb_rlogMatrix_Unique_topVarGenes_rownames[,3] + All_clusters_25kb_rlogMatrix_Unique_topVarGenes_rownames[,8])/2,
           (All_clusters_25kb_rlogMatrix_Unique_topVarGenes_rownames[,4] + All_clusters_25kb_rlogMatrix_Unique_topVarGenes_rownames[,9])/2,
           (All_clusters_25kb_rlogMatrix_Unique_topVarGenes_rownames[,5] + All_clusters_25kb_rlogMatrix_Unique_topVarGenes_rownames[,10])/2,
           (All_clusters_25kb_rlogMatrix_Unique_topVarGenes_rownames[,6] + All_clusters_25kb_rlogMatrix_Unique_topVarGenes_rownames[,11])/2,
           All_clusters_25kb_rlogMatrix_Unique_topVarGenes_rownames[,15])
colnames(All_clusters_25kb_rlogMatrix_Unique_topVarGenes_Avg) <- c("ENSGID", "Avg_0", "Avg_30", "Avg_60", "Avg_120", "Avg_24", "Cluster") 

# Add StDev column
All_clusters_25kb_rlogMatrix_Unique_topVarGenes_Avg$StDev <- apply(All_clusters_25kb_rlogMatrix_Unique_topVarGenes_Avg[,2:6], 1, sd)

# Center and standardize 

All_clusters_25kb_rlogMatrix_Unique_topVarGenes_Avg_Standardized <- data.frame(All_clusters_25kb_rlogMatrix_Unique_topVarGenes_Avg[,1], 
           (All_clusters_25kb_rlogMatrix_Unique_topVarGenes_Avg[,2:6] - rowMeans(All_clusters_25kb_rlogMatrix_Unique_topVarGenes_Avg[,2:6]))/All_clusters_25kb_rlogMatrix_Unique_topVarGenes_Avg$StDev,
           All_clusters_25kb_rlogMatrix_Unique_topVarGenes_Avg[,7])
colnames(All_clusters_25kb_rlogMatrix_Unique_topVarGenes_Avg_Standardized) <- c("ENSEMBL", "0", "0.5", "1", "2", "24", "Cluster")

# Change to long format
gather(All_clusters_25kb_rlogMatrix_Unique_topVarGenes_Avg_Standardized, "Time", "RNA_reads", 2:6) -> All_clusters_25kb_rlogMatrix_Unique_topVarGenes_Avg_Standardized_Tidy          

######### Create table of meth, accessibility and RNA observations together ######### 

RNA_table <- All_clusters_25kb_rlogMatrix_Unique_topVarGenes_Avg_Standardized_Tidy[,2:4]
colnames(RNA_table) <- c("Cluster", "Time", "Value")
RNA_table$Type <- "DEseq_counts"
ATAC_table <- TCseq_Clustering_Counts_Table_Standardized_tidy[,5:7]
colnames(ATAC_table) <- c("Cluster", "Time", "Value")
ATAC_table$Type <- "ATAC_counts"
Meth_table <- meth_time_table_centered_longformat_filter5[,7:9]
colnames(Meth_table) <- c("Cluster", "Time", "Value")
Meth_table$Type <- "METH"

# Bind tables together
combined_table <- rbind(RNA_table, ATAC_table, Meth_table)
combined_table$Time <- as.factor(combined_table$Time)

# Make plot of all three box plot series
Cluster_combined_boxplot <- ggplot(filter(combined_table, Cluster %in% c(5)), aes(x=Time, y=Value)) +
  geom_boxplot(outlier.shape = NA, aes(fill=Type), position=position_dodge(0.60), width=0.4)+
  stat_summary(fun.y=mean, geom="smooth", aes(group=Type, color=Type))+
  xlab("Time (hrs)")+
  ylab("Standard Difference")+
  ggtitle("Cluster 5 filter5_CpG")+
  theme_bw()
Cluster_combined_boxplot

# Cluster_access_boxplot + coord_cartesian(ylim = c(-3,3))

### Confirm if CpGs should be depth filtered when calculating average methlyation for the region

# Read in data
Depth1_0_avg_table <- read.delim2(file = "cluster_allCluster_0_depth1.roiavg.table", header = FALSE)
Depth5_0_avg_table <- read.delim2(file = "cluster_allCluster_0_depth5.roiavg.table", header = FALSE)
Depth1_120_avg_table <- read.delim2(file = "cluster_allCluster_120_depth1.roiavg.table", header = FALSE)
Depth5_120_avg_table <- read.delim2(file = "cluster_allCluster_120_depth5.roiavg.table", header = FALSE)

# Change data type of methylation columns
Depth1_0_avg_table$V5 <- as.numeric(as.character(Depth1_0_avg_table$V5))
Depth5_0_avg_table$V5 <- as.numeric(as.character(Depth5_0_avg_table$V5))
Depth1_120_avg_table$V5 <- as.numeric(as.character(Depth1_120_avg_table$V5))
Depth5_120_avg_table$V5 <- as.numeric(as.character(Depth5_120_avg_table$V5))

# Merge tables

DepthCompare_0_avg_table <- merge(Depth1_0_avg_table, Depth5_0_avg_table, by = c("V1", "V2", "V3"))
DepthCompare_120_avg_table <- merge(Depth1_120_avg_table, Depth5_120_avg_table, by = c("V1", "V2", "V3"))

# Plot with ggplot 

ggplot(DepthCompare_0_avg_table, aes(x = DepthCompare_0_avg_table$V5.x, y = DepthCompare_0_avg_table$V5.y))+
  geom_point(alpha=0.15, color="firebrick4")+
  ggtitle("Read Depth Comparison for Meth Avg Calculation 0hr")+
  ylab("Avg Region Methylation Depth>=5")+
  xlab("Avg Region Methylation Depth>=1")+
  ylim(0, 1)+
  xlim(0, 1)+
  theme(plot.title = element_text(lineheight=.8, face="bold", vjust=1),
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust=0.3))+ theme_bw()

ggplot(DepthCompare_120_avg_table, aes(x = DepthCompare_120_avg_table$V5.x, y = DepthCompare_120_avg_table$V5.y))+
  geom_point(alpha=0.15, color="firebrick4")+
  ggtitle("Read Depth Comparison for Meth Avg Calculation 2hr")+
  ylab("Avg Region Methylation Depth>=5")+
  xlab("Avg Region Methylation Depth>=1")+
  ylim(0, 1)+
  xlim(0, 1)+
  theme(plot.title = element_text(lineheight=.8, face="bold", vjust=1),
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust=0.3))+ theme_bw()


### Create table of average methylation across time 

# Load required library
library(readr)

# Read in data and add column names
TimePoint_0_avg_table <- read_tsv(file = "cluster_allCluster_0_depth1.roiavg.table", col_names = FALSE)[,c(1,2,3,4,9)]
colnames(TimePoint_0_avg_table) <- c("Chr", "Start", "End", "PeakID", "Time0")

TimePoint_30_avg_table <- read_tsv(file = "cluster_allCluster_30_depth1.roiavg.table", col_names = FALSE)[,c(1,2,3,4,9)]
colnames(TimePoint_30_avg_table) <- c("Chr", "Start", "End", "PeakID", "Time0.5")

TimePoint_60_avg_table <- read_tsv(file = "cluster_allCluster_60_depth1.roiavg.table", col_names = FALSE)[,c(1,2,3,4,9)]
colnames(TimePoint_60_avg_table) <- c("Chr", "Start", "End", "PeakID", "Time1")

TimePoint_120_avg_table <- read_tsv(file = "cluster_allCluster_120_depth1.roiavg.table", col_names = FALSE)[,c(1,2,3,4,9)]
colnames(TimePoint_120_avg_table) <- c("Chr", "Start", "End", "PeakID", "Time2")

TimePoint_24hr_avg_table <- read_tsv(file = "cluster_allCluster_24hr_depth1.roiavg.table", col_names = FALSE)[,c(1,2,3,4,9)]
colnames(TimePoint_24hr_avg_table) <- c("Chr", "Start", "End", "PeakID", "Time24")

# Make a list of all data frame
TimeCourse_AvgMeth_List <- list(TimePoint_0_avg_table, TimePoint_30_avg_table, TimePoint_60_avg_table, TimePoint_120_avg_table, TimePoint_24hr_avg_table)

# Use reduce to iteratively merge all data frames
# Results in only getting loci for which an average methylation value was able to be calculated at ALL time points
# Regions that have no CpG or no covered CpG will not be in this table 
TimeCourse_AvgMeth_Table <- Reduce(function(x, y) merge(x, y, all=F, 
                                    by=c("Chr", "Start", "End", "PeakID")), TimeCourse_AvgMeth_List, accumulate=F)

# Same table as above (TimeCourse_AvgMeth_Table) but now contains NAs if there was no data at a particular time point instead of excluding that region
TimeCourse_AvgMeth_Table_NAs <- Reduce(function(x, y) merge(x, y, all=T, 
                                                        by=c("Chr", "Start", "End", "PeakID")), TimeCourse_AvgMeth_List, accumulate=F)
# Join this table with RNA and ATAC tables 
