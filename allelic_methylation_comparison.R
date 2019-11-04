# Load required libraries
library(dplyr)
library(ggplot2)
library(gghighlight)
library(Cairo)

### Compare CpG methylation between ATACMe and WGBS within verified allelicly methylated regions

# Read in bed file of CpGs
ATACMe_VerifiedAMR_CpGs <- read.delim2("ATACMe_VerifiedAMR_CpGs.bed", header = FALSE)
WGBS0hr_VerifiedAMR_CpGs <- read.delim2("WGBS0hr_VerifiedAMR_CpGs.bed", header = FALSE)
WGBS24hr_VerifiedAMR_CpGs <- read.delim2("WGBS24hr_VerifiedAMR_CpGs.bed", header = FALSE)

# Convert to tbl 
ATACMe_VerifiedAMR_CpGs <- tbl_df(ATACMe_VerifiedAMR_CpGs)
WGBS0hr_VerifiedAMR_CpGs <- tbl_df(WGBS0hr_VerifiedAMR_CpGs)
WGBS24hr_VerifiedAMR_CpGs <- tbl_df(WGBS24hr_VerifiedAMR_CpGs)

# Convert meth column from factor to numeric 
ATACMe_VerifiedAMR_CpGs$V4 <- as.numeric(as.character(ATACMe_VerifiedAMR_CpGs$V4))
WGBS0hr_VerifiedAMR_CpGs$V4 <- as.numeric(as.character(WGBS0hr_VerifiedAMR_CpGs$V4))
WGBS24hr_VerifiedAMR_CpGs$V4 <- as.numeric(as.character(WGBS24hr_VerifiedAMR_CpGs$V4))

# Filter for CpGs that have at least 5 reads 
ATACMe_VerifiedAMR_CpGs_ReadFilter5 <- filter(ATACMe_VerifiedAMR_CpGs, V5 >= 5)
WGBS0hr_VerifiedAMR_CpGs_ReadFilter5 <- filter(WGBS0hr_VerifiedAMR_CpGs, V5 >= 5)
WGBS24hr_VerifiedAMR_CpGs_ReadFilter5 <- filter(WGBS24hr_VerifiedAMR_CpGs, V5 >= 5)

# Only use CpGs that have at least 5 reads in all three libraries 
ATACMe_and_WGBS0_SharedCpGs_ReadDepth5 <- merge(ATACMe_VerifiedAMR_CpGs_ReadFilter5, WGBS0hr_VerifiedAMR_CpGs_ReadFilter5, by = c("V1", "V2", "V3"), all = FALSE)
ATACMe_and_WGBS24_SharedCpGs_ReadDepth5 <- merge(ATACMe_VerifiedAMR_CpGs_ReadFilter5, WGBS24hr_VerifiedAMR_CpGs_ReadFilter5, by = c("V1", "V2", "V3"), all = FALSE)

# Filter for CpGs that actually show intermediate methylation (0.40 to 0.60) in WGBS dataset
ATACMe_and_WGBS0_SharedCpGs_ReadDepth5_intermediateMeth <- filter(ATACMe_and_WGBS0_SharedCpGs_ReadDepth5, V4.y > 0.35 & V4.y < 0.65)
ATACMe_and_WGBS24_SharedCpGs_ReadDepth5_intermediateMeth <- filter(ATACMe_and_WGBS24_SharedCpGs_ReadDepth5, V4.y > 0.35 & V4.y < 0.65)

# Correlation Test of intermediate methylation CpGs (as calculated in WGBS dataset)
cor.test(ATACMe_and_WGBS0_SharedCpGs_ReadDepth5_intermediateMeth$V4.x, ATACMe_and_WGBS0_SharedCpGs_ReadDepth5_intermediateMeth$V4.y, method=c("pearson"))
cor.test(ATACMe_and_WGBS24_SharedCpGs_ReadDepth5_intermediateMeth$V4.x, ATACMe_and_WGBS24_SharedCpGs_ReadDepth5_intermediateMeth$V4.y, method=c("pearson"))

# Correlation test of any CpG that fell within an AMR and had sufficient coverage in both ATACMe and WGBS
cor.test(ATACMe_and_WGBS0_SharedCpGs_ReadDepth5$V4.x, ATACMe_and_WGBS0_SharedCpGs_ReadDepth5$V4.y, method=c("pearson"))
cor.test(ATACMe_and_WGBS24_SharedCpGs_ReadDepth5$V4.x, ATACMe_and_WGBS24_SharedCpGs_ReadDepth5$V4.y, method=c("pearson"))

# Scatterplots of any CpG that fell within an AMR and had sufficient coverage in both ATACMe and WGBS
ATACMe_WGBS_AMR_CpG_Scatterplot_0 <- ggplot(ATACMe_and_WGBS0_SharedCpGs_ReadDepth5, aes(x=V4.x, y=V4.y)) +
  geom_point(size=4.5, shape=20, alpha=.35)+
  xlab("ATACMe CpG Methylation Level")+
  ylab("WGBS CpG Methylation Level (0hr)")+
  ggtitle("CpG Methylation Levels within Validated AMRs")+
  theme_bw()

ggsave(filename = "CpG_Methylation_Scatterplot_Validated_AMRs_0hr.eps",
       plot = print(ATACMe_WGBS_AMR_CpG_Scatterplot_0),
       device = cairo_ps)

ATACMe_WGBS_AMR_CpG_Scatterplot_24 <- ggplot(ATACMe_and_WGBS24_SharedCpGs_ReadDepth5, aes(x=V4.x, y=V4.y)) +
  geom_point(size=4.5, shape=20, alpha=.35)+
  xlab("ATACMe CpG Methylation Level")+
  ylab("WGBS CpG Methylation Level (24hr)")+
  ggtitle("CpG Methylation Levels within Validated AMRs")+
  theme_bw()

ggsave(filename = "CpG_Methylation_Scatterplot_Validated_AMRs_24hr.eps",
       plot = print(ATACMe_WGBS_AMR_CpG_Scatterplot_24),
       device = cairo_ps)


### Compare CpG methylation between ATACMe and WGBS within peak regions

# Read in bed file of CpGs
ATACMe_gfTCseqPeaks_CpGs <- read.delim2("ATACMe_gfTCseqPeaks_CpGs.bed", header = FALSE)
WGBS0hr_gfTCseqPeaks_CpGs <- read.delim2("WGBS0hr_gfTCseqPeaks_CpGs.bed", header = FALSE)
WGBS24hr_gfTCseqPeaks_CpGs <- read.delim2("WGBS24hr_gfTCseqPeaks_CpGs.bed", header = FALSE)

# Convert to tbl 
ATACMe_gfTCseqPeaks_CpGs <- tbl_df(ATACMe_gfTCseqPeaks_CpGs)
WGBS0hr_gfTCseqPeaks_CpGs <- tbl_df(WGBS0hr_gfTCseqPeaks_CpGs)
WGBS24hr_gfTCseqPeaks_CpGs <- tbl_df(WGBS24hr_gfTCseqPeaks_CpGs)

# Convert meth column from factor to numeric 
ATACMe_gfTCseqPeaks_CpGs$V4 <- as.numeric(as.character(ATACMe_gfTCseqPeaks_CpGs$V4))
WGBS0hr_gfTCseqPeaks_CpGs$V4 <- as.numeric(as.character(WGBS0hr_gfTCseqPeaks_CpGs$V4))
WGBS24hr_gfTCseqPeaks_CpGs$V4 <- as.numeric(as.character(WGBS24hr_gfTCseqPeaks_CpGs$V4))

# Filter for CpGs that have at least 5 reads 
ATACMe_gfTCseqPeaks_CpGs_ReadFilter5 <- filter(ATACMe_gfTCseqPeaks_CpGs, V5 >= 5)
WGBS0hr_gfTCseqPeaks_CpGs_ReadFilter5 <- filter(WGBS0hr_gfTCseqPeaks_CpGs, V5 >= 5)
WGBS24hr_gfTCseqPeaks_CpGs_ReadFilter5 <- filter(WGBS24hr_gfTCseqPeaks_CpGs, V5 >= 5)

# Only use CpGs that have at least 5 reads between ATACMe and WGBS 
ATACMe_and_WGBS0_SharedCpGs_gfTCseqPeaks_ReadDepth5 <- merge(ATACMe_gfTCseqPeaks_CpGs_ReadFilter5, WGBS0hr_gfTCseqPeaks_CpGs_ReadFilter5, by = c("V1", "V2", "V3"), all = FALSE)
ATACMe_and_WGBS24hr_SharedCpGs_gfTCseqPeaks_ReadDepth5 <- merge(ATACMe_gfTCseqPeaks_CpGs_ReadFilter5, WGBS24hr_gfTCseqPeaks_CpGs_ReadFilter5, by = c("V1", "V2", "V3"), all = FALSE)

# Randomly sample lines to make plotting easier 
ATACMe_and_WGBS0_SharedCpGs_gfTCseqPeaks_ReadDepth5_downSample <- sample_n(ATACMe_and_WGBS0_SharedCpGs_gfTCseqPeaks_ReadDepth5, 1800, replace = FALSE, weight = NULL, .env = NULL)
ATACMe_and_WGBS24hr_SharedCpGs_gfTCseqPeaks_ReadDepth5_downSample <- sample_n(ATACMe_and_WGBS24hr_SharedCpGs_gfTCseqPeaks_ReadDepth5, 1800, replace = FALSE, weight = NULL, .env = NULL)

# Correlation Test 
cor.test(ATACMe_and_WGBS0_SharedCpGs_gfTCseqPeaks_ReadDepth5_downSample$V4.x, ATACMe_and_WGBS0_SharedCpGs_gfTCseqPeaks_ReadDepth5_downSample$V4.y, method=c("pearson"))
cor.test(ATACMe_and_WGBS24hr_SharedCpGs_gfTCseqPeaks_ReadDepth5_downSample$V4.x, ATACMe_and_WGBS24hr_SharedCpGs_gfTCseqPeaks_ReadDepth5_downSample$V4.y, method=c("pearson"))


# Scatterplot comparison 
ATACMe_WGBS_Peaks_CpG_Scatterplot_0 <- ggplot(ATACMe_and_WGBS0_SharedCpGs_gfTCseqPeaks_ReadDepth5_downSample, aes(x=V4.x, y=V4.y)) +
  geom_point(size=4.5, shape=20, alpha=.35)+
  xlab("ATACMe CpG Methylation Level")+
  ylab("WGBS CpG Methylation Level (0hr)")+
  ggtitle("CpG Methylation Levels within Peak Regions")+
  theme_bw()

ggsave(filename = "CpG_Methylation_Scatterplot_PeakRegions_0hr.eps",
       plot = print(ATACMe_WGBS_Peaks_CpG_Scatterplot_0),
       device = cairo_ps)

ATACMe_WGBS_Peaks_CpG_Scatterplot_24 <- ggplot(ATACMe_and_WGBS24hr_SharedCpGs_gfTCseqPeaks_ReadDepth5_downSample, aes(x=V4.x, y=V4.y)) +
  geom_point(size=4.5, shape=20, alpha=.35)+
  xlab("ATACMe CpG Methylation Level")+
  ylab("WGBS CpG Methylation Level (24hr)")+
  ggtitle("CpG Methylation Levels within Peak Regions")+
  theme_bw()

ggsave(filename = "CpG_Methylation_Scatterplot_PeakRegions_24hr.eps",
       plot = print(ATACMe_WGBS_Peaks_CpG_Scatterplot_24),
       device = cairo_ps)






