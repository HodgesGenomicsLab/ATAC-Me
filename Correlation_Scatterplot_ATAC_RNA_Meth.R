library(dplyr)
library(ggplot2)
library(tidyr)


#### Create scatterplot of accessibility versus associated gene expression #### 

# Use "TCseq_Clustering_Counts_Table" from ATAC_and_Meth_Table_Construction.RData
# Merge above table with RNA data contained in "ClusterAnnotation_25kb_withCounts" from ATAC_Cluster_GeneExpression_Analysis_EH4889.RData
# User normalized counts for both ATAC and RNA data

ATAC_and_RNA_Counts <- merge(TCseq_Clustering_Counts_Table, ClusterAnnotation_25kb_withCounts, by.x = "PeakID", by.y = "V4", all = FALSE)
ATAC_and_RNA_Counts2 <- merge(All_clusters_25kb_rlogMatrix_Unique_topVarGenes_Avg, ATAC_and_RNA_Counts, by.x = "ENSGID", by.y = "ENSEMBL", all = FALSE)

ATAC_and_RNA_Counts2[ ,1:18]
ATAC_and_RNA_Counts2_unique <- unique(ATAC_and_RNA_Counts2[ ,1:18])

# Change column names 
colnames(ATAC_and_RNA_Counts2_unique) <- c("ENSGID", "RNA_0", "RNA_30", "RNA_60", "RNA_120","RNA_24hr", "Cluster1", "StDev", "PeakID", "ATAC_0", "ATAC_30", "ATAC_60", "ATAC_120", "ATAC_24", "Chr", "Start", "End", "Cluster2")

# Melt data 
RNA_melted_data <- gather(ATAC_and_RNA_Counts2_unique[ ,1:6], "RNA_Time", "RNA Counts", 2:6) 
ATAC_melted_data <- gather(ATAC_and_RNA_Counts2_unique[ ,9:18], "ATAC_Time", "ATAC Counts", 2:6)

# Make new data frame of melted values
RNA_and_ATAC_melted_data <- cbind(RNA_melted_data$`RNA Counts`, ATAC_melted_data$`ATAC Counts`)
colnames(RNA_and_ATAC_melted_data) <- c("RNA", "ATAC")
RNA_and_ATAC_melted_data_tibble <- as_tibble(RNA_and_ATAC_melted_data)

RNA_and_ATAC_melted_data_tibble <- mutate(RNA_and_ATAC_melted_data_tibble, ATAC_log2 = log2(ATAC))

# Correlation test of data
cor.test(RNA_and_ATAC_melted_data_tibble$RNA, RNA_and_ATAC_melted_data_tibble$ATAC_log2, method=c("pearson"))


# Make scatterplot
ATAC_and_RNA_Counts2_unique_plot <- ggplot(RNA_and_ATAC_melted_data_tibble, aes(x=ATAC_log2, y=RNA))+
  ggtitle("RNA vs. Access")+
  geom_point(alpha=.80, color="indianred3")+
  xlab("log2(ATAC Read Counts)")+
  ylab("rlog(RNA Read Counts)")+
  theme(plot.title = element_text(lineheight=.8, face="bold", vjust=1),
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust=0.3))+ theme_bw()

# Save scatterplot as eps 
# Use ggsave command below, otherwise transparency will not be maintained
library(Cairo)
ggsave(filename = "ATAC_versus_RNA_Correlation_Scatterplot.eps",
       plot = print(ATAC_and_RNA_Counts2_unique_plot),
       device = cairo_ps)

#### Remake scatterplot using log2 transformed normalized counts rather than rlog transformed values

# Pull out counts and add column that is ensembl ID without decimal identifier
counts_ddsTC2 <- data.frame(ENSEMBL = gsub("\\..*","", rownames(counts_ddsTC)), counts_ddsTC)

# Get average of replicates 
counts_ddsTC2_avg <- data.frame(ENSGID = counts_ddsTC2[,1], avg_0 = ((counts_ddsTC2[,2] + counts_ddsTC2[,7])/2),
           avg_30 = ((counts_ddsTC2[,3] + counts_ddsTC2[,8])/2),
           avg_60 = ((counts_ddsTC2[,4] + counts_ddsTC2[,9])/2),
           avg_120 = ((counts_ddsTC2[,5] + counts_ddsTC2[,10])/2),
           avg_24 = ((counts_ddsTC2[,6] + counts_ddsTC2[,11])/2))
# Merge previous table of ATAC and RNA counts 

ATAC_and_RNA_Counts2_unique_withNormCounts <- merge(ATAC_and_RNA_Counts2_unique, counts_ddsTC2_avg, by = "ENSGID", all = FALSE)

# Melt data 
RNA_melted_data2 <- gather(ATAC_and_RNA_Counts2_unique_withNormCounts[ ,c(1,19:23)], "RNA_Time", "RNA Counts", 2:6) 
ATAC_melted_data2 <- gather(ATAC_and_RNA_Counts2_unique_withNormCounts[ ,9:18], "ATAC_Time", "ATAC Counts", 2:6)

# Make new data frame of melted values
RNA_and_ATAC_melted_data2 <- cbind(RNA_melted_data2$`RNA Counts`, ATAC_melted_data2$`ATAC Counts`)
colnames(RNA_and_ATAC_melted_data2) <- c("RNA", "ATAC")
RNA_and_ATAC_melted_data_tibble2 <- as_tibble(RNA_and_ATAC_melted_data2)

# Mutate to transform the data values into log2 form
RNA_and_ATAC_melted_data_tibble2 <- mutate(RNA_and_ATAC_melted_data_tibble2, ATAC_log2 = log2(ATAC), RNA_log2 = log2(RNA))

#Scatterplot
ATAC_and_RNA_Counts2_unique_plot2 <- ggplot(RNA_and_ATAC_melted_data_tibble2, aes(x=ATAC_log2, y=RNA_log2))+
  ggtitle("RNA vs. Access")+
  geom_point(alpha=0.40, color="steelblue")+
  xlab("log2(ATAC Read Counts)")+
  ylab("log2(RNA Read Counts)")+
  theme(plot.title = element_text(lineheight=.8, face="bold", vjust=1),
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust=0.3))+ theme_bw()



#### Create scatterplot of DNA Methylation versus associated gene expression #### 
# Get average DNA methylation of the regions rather than use CpGs 

# Load required library
library(readr)

# Set directory where files read in are
setwd("/Users/kellybarnett/Desktop/EH4889_noshift")

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

# Merge meth table with chromatin accessibility and RNA data

ATAC_RNA_Meth_Table <- merge(ATAC_and_RNA_Counts2_unique, TimeCourse_AvgMeth_Table, by = c("Chr", "Start", "End", "PeakID"), all=F)

# Melt data 
RNA_melted_data2 <- gather(ATAC_RNA_Meth_Table[ ,1:11], "RNA_Time", "RNA Counts", 6:10) 
Meth_melted_data <- gather(ATAC_RNA_Meth_Table[ ,c(1:5, 11, 19:23)], "Meth_Time", "Avg Meth", 7:11)

# Make new data frame of melted values
RNA_and_METH_melted_data <- cbind(RNA_melted_data2$`RNA Counts`, Meth_melted_data$`Avg Meth`)
colnames(RNA_and_METH_melted_data) <- c("RNA", "METH")
RNA_and_METH_melted_data <- as_tibble(RNA_and_METH_melted_data)

# Correlation test of data
cor.test(RNA_and_METH_melted_data$RNA, RNA_and_METH_melted_data$METH, method=c("pearson"))

ATAC_and_METH_plot <- ggplot(RNA_and_METH_melted_data, aes(x=METH, y=RNA))+
  ggtitle("RNA vs. Meth")+
  geom_point(alpha=0.80, color="royalblue3")+
  xlab("Avg DNA Meth for Region")+
  ylab("rlog transformed RNA read counts")+
  theme(plot.title = element_text(lineheight=.8, face="bold", vjust=1),
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust=0.3))+ theme_bw()

# Save plot as eps
library(Cairo)
ggsave(filename = "METH_versus_RNA_Correlation_Scatterplot.eps",
       plot = print(ATAC_and_METH_plot),
       device = cairo_ps)

