# Peaks directory
# /data/hodges_lab/EH4889/trimmed_reads/mapped_reads/methprocessed/macs/broadpeaks

# Peak files
# Files have been filtered for q-value
# Only includes peaks reproducible across replicates
# THP1_ATAC-Me_0_concordant.broadPeak.q10
# THP1_ATAC-Me_30_concordant.broadPeak.q10
# THP1_ATAC-Me_60_concordant.broadPeak.q10
# THP1_ATAC-Me_120_concordant.broadPeak.q10
# THP1_ATAC-Me_24hr_concordant.broadPeak.q10

# Bam file directory
# /data/hodges_lab/EH4889/trimmed_reads/mapped_reads/methprocessed/macs

# Bam file list 
# R1_ATAC-Me_THP1_PMA_Rep1A_0.hum		R1_ATAC-Me_THP1_PMA_Rep2A_0.hum
# R1_ATAC-Me_THP1_PMA_Rep1A_30.hum	R1_ATAC-Me_THP1_PMA_Rep2A_30.hum
# R1_ATAC-Me_THP1_PMA_Rep1A_60.hum	R1_ATAC-Me_THP1_PMA_Rep1A_60.hum
# R1_ATAC-Me_THP1_PMA_Rep1A_120.hum	R1_ATAC-Me_THP1_PMA_Rep1A_120.hum
# R1_ATAC-Me_THP1_PMA_Rep1A_24hr.hum	R1_ATAC-Me_THP1_PMA_Rep1A_24hr.hum

# Set path to R-package libraries
.libPaths("/data/hodges_lab/R_3.6.0_VizPortal/") 
# Load TCseq package
library(TCseq)

# Set directory location of peak files
dir <- "/data/hodges_lab/EH4889/trimmed_reads/mapped_reads/methprocessed/macs/broadpeaks"

# This defines the common genomic intervals for analysis across the entire time course
# TCseq merges overlapping regions into a contiguous genomic interval
# use the pattern argument to match the files of interest in the directory
gf <- peakreference(dir = dir, pattern = "concordant.broadPeak.q10") 

# Set vector of unique sample ids
sampleid = c("s1", "s2", "s3", "s4","s5", "s6", "s7", "s8", "s9","s10")

# Set vector defining the time point conditions
timepoint = c("0min","30min","60min","120min","24h")

# Set vector defining the experimental treatments
group = c(1, 2, 3, 4, 5)

# Vector list of the bam files
BAMfile =c("R1_ATAC-Me_THP1_PMA_Rep1A_0.hum","R1_ATAC-Me_THP1_PMA_Rep1A_30.hum","R1_ATAC-Me_THP1_PMA_Rep1A_60.hum","R1_ATAC-Me_THP1_PMA_Rep1A_120.hum","R1_ATAC-Me_THP1_PMA_Rep1A_24hr.hum","R1_ATAC-Me_THP1_PMA_Rep2A_0.hum","R1_ATAC-Me_THP1_PMA_Rep2A_30.hum","R1_ATAC-Me_THP1_PMA_Rep2A_60.hum","R1_ATAC-Me_THP1_PMA_Rep2A_120.hum","R1_ATAC-Me_THP1_PMA_Rep2A_24hr.hum") 

# Build data.frame of experimental design
exp_bamfile = data.frame(sampleid,timepoint,group,BAMfile)

# Input genomic regions and experimental design into TCA object
tca <- TCA(design = exp_bamfile, genomicFeature = gf)

# Set directory location of bam files
dir.BAM <- "/data/hodges_lab/EH4889/trimmed_reads/mapped_reads/methprocessed/macs"

# Get count of reads overlapping peak regions of interest
tca <- countReads(tca, dir = dir.BAM)

# Run analysis for differential events across the time course
# Uses GLM based analysis from edgeR
# Filters out regions that have low count values
tca <- DBanalysis(tca, categories = "timepoint", norm.lib = TRUE, filter.type = "cpm", filter.value = 5, samplePassfilter = 2)

# Builds table of differential time course events that pass the defined pvalue filtering
tca <- timecourseTable(tca, value = "expression", norm.method = "rpkm", filter = TRUE, pvalue = "fdr", pvalue.threshold = 0.005, direction = "both")

# Write table of values for differential time course events
t <- tcTable(tca)
write.table(t, sep="\t", file="TCseq_EH4889_ATAC_ExpressionTable.txt")
save.image(file = "TCseq_EH4889_ATAC_noShift.RData)

# Cluster time course data into 5 clusters via fuzzy c-means clustering
tca <- timeclust(tca, algo = "cm", k = 5, standardize = TRUE)

# Make plot of timecourse clustering 
p <- timeclustplot(tca, value = "z-score(RPKM)", cols = 3, membership.color = topo.colors(30))

# Print out plot to pdf
pdf(file = "TCPlot_q10_broad_concord_noshift_cpm5.pdf", width = 10, height = 7)

p

dev.off()

## Code to output bed files of peaks (created by the TCSeq
## peakreference() function) belonging to each cluster

## Peaks whose accessibility did not significantly change between any
## timepoints are not included in the output

## Make a dataframe of peaks with cluster information
# at this point, tca is the object that resulted from calling timeclust()
cluster = tca@clusterRes@cluster
tca_clust_gf = tca@genomicFeature
rownames(tca_clust_gf) = tca_clust_gf[,4]
gf_clust = tca_clust_gf[names(cluster),]
gf_clust = cbind(gf_clust,cluster)
# gf_clust is a dataframe whose columns include chr, start, end, the
Peak ID given by TCSeq, and the cluster assigned by the clustering
algorithm of choice
## Write bed files per cluster
for(i in unique(gf_clust$cluster)){
write.table(gf_clust[gf_clust$cluster==i,],paste("peak_subset_cluster_",i,".bed",sep=""),sep="\t",row.names=FALSE,quote=FALSE,col.names=FALSE)
}
## The resulting bed file needs to be sorted


# Write gf variable to bed file. This is a list of all input peak regions after TC_seq runs a union merge
write.table(gf, file="gf_tcseq.bed", quote=F, sep="\t", row.names=F, col.names=F)
