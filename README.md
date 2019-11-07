# ATAC-Me
Scripts and analysis approaches for ATAC-Me origin data
Most scripts assume access to a high power computer or cluster computing environment

## Software 

    - bedtools v2.26.0
    - samtools 1.5
    - deeptools 3.1.2
    - STAR_2.6.1a 
    - WALT v1.0 
    - bowtie2 version 2.2.6 
    - FastQC v0.11.4
    - cutadapt 1.8.3 
    - trim_galore 0.4.0 
    - DeSeq2 1.18.1 (R version 3.4.3)
    - ChIPeeker 1.14.2 (R version 3.4.3)
    - clusterProfiler 3.6.0 (R version 3.4.3)
    - bagfoot 0.9.6 (R version 3.4.3)

## General Processing Scripts

1. Trimming
    - trimming_loop.sh
    - trim.slrm
2. Mapping 
    - mapping_loop_walt.sh
    - mapping_walt.slrm
3. Methylation levels
    - methpipe_loop.sh
    - methpipe.slrm
4. Peak calling
    - macs2_loop.sh
    - mapping_loop_walt.sh
    - peak_filtering.txt
    
## Figure 1
###### ATAC-Me vs. Standard ATAC Read Count Correlation
    - read_count_correlation.txt
###### Allelic Methylation Analysis
    - allelic_methylation_comparison.R

## Figure 2

###### TC-seq Processing and Clustering
    - TCseq_Analysis.R

## Figure 3

###### DNA methylation Heatmaps Across TCseq Cluster Regions
    - TCseq_cluster_meth_heatmaps.txt

## Figure 4
###### KEGG Pathway/GO Analysis of Dynamic ATAC Peak Neighboring Genes

###### K-means Clustering of RNA-seq Differential Genes
    - DEseq2_RNA_Timecourse.R
###### Hierarchical Clustering of Top Variable Genes Near Dynamic ATAC Peaks + Methylation
    - topVariable_Genes_Heatmap.R 
## Figure 5
###### Scatterplot Correlations of ATAC vs. RNA and Methylation vs. RNA
    - Correlation_Scatterplot_ATAC_RNA_Meth.R
###### Boxplot Time Series Integrating ATAC, RNA, Methylation
    - Integrative_Boxplot_TimeSeries.R
## Figure 6
###### Boxplot Time Series for Extended Time Point Data
    - Extended_Timepoints_Boxplot_Comparision.R
