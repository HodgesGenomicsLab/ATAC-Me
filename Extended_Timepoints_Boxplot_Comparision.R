library(tidyverse)
library(Cairo)
library(ggpubr)

setwd("~/Desktop/3854_Vantage/Meth_bed_files/")

# Build data table for dynamic region CpGs
files_dyn <- list.files(pattern = ".*DynamicOnly.*bed")
df_list_dyn <- lapply(files_dyn, read_tsv, col_names = c("Chr", "Start", "End", "Meth", "Reads", "Chr_TC", "Start_TC", "End_TC", "PeakID", "Cluster", "Overlap"))
names(df_list_dyn) <- c("0hr","24hr", "48hr", "72hr", "72hr+VC")
df_list_dyn_tidy <- bind_rows(df_list_dyn, .id = "Sample") %>% select(-Overlap) %>% mutate_at("Sample", as.factor)

# Build data table for static regions CpGs
files_stat <- list.files(pattern = ".*StaticOnly.*bed")
df_list_stat <- lapply(files_stat, read_tsv, col_names = c("Chr", "Start", "End", "Meth", "Reads", "Chr_TC", "Start_TC", "End_TC", "PeakID", "Overlap", "Cluster"))
names(df_list_stat) <- c("0hr","24hr", "48hr", "72hr", "72hr+VC")
df_list_stat_tidy <- bind_rows(df_list_stat, .id = "Sample") %>% select(-Overlap) %>% mutate_at("Sample", as.factor)

# Build ggplot boxplot for dynamic regions 
Boxplot_Meth_OpeningRegions <- df_list_dyn_tidy %>% 
  filter(Cluster != c(2,1)) %>% 
  ggplot(aes(x = Sample, y = Meth)) + 
  geom_boxplot(color = "black", fill = "grey70", width = 0.50,
               outlier.alpha = 0.01, outlier.size = 1.5, outlier.shape = 1, outlier.colour = "firebrick4", outlier.stroke = 0.5) + 
  ggtitle("CpG Methylation Among Chromatin Opening Regions") +
  ylab("CpG Methylation Fraction") +
  xlab("Treatment Group (hrs of PMA exposure)") +
  theme_bw() +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=14)) +  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif")

# Get average difference between 72hr and 72hr+VitC

df_list_dyn_tidy %>% group_by(Sample) %>% 
  summarize(meth_mean = mean(Meth)) %>% 
  head()






# Build plot with ggpubr
my_comparisons <- list( c("0hr", "24hr"), c("24hr", "48hr"), c("24hr", "72hr"), c("72hr", "72hr+VC"))

df_list_dyn_tidy %>% 
  filter(Cluster != c(2,1)) %>% 
  ggboxplot(x = "Sample", y = "Meth", color = "black", fill = "grey70", width = 0.50,
            ylab = "CpG Methylation Fraction", xlab = "Treatment Group (hrs of PMA exposure)", title = "CpG Methylation Among Regions Open at 24hrs", outlier.alpha = 0.01) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif")
 


Boxplot_Meth_OpeningRegions + compare_means(data = df_list_dyn_tidy, method = "wilcox.test", paired = FALSE)

# Test pairs of regions for significance 
wilcox.test(filter(df_list_dyn_tidy, Sample == "72hr")$Meth, filter(df_list_dyn_tidy, Sample == "72hr+VC")$Meth, alternative = "two.sided")
wilcox.test(filter(df_list_dyn_tidy, Sample == "24hr")$Meth, filter(df_list_dyn_tidy, Sample == "48hr")$Meth, alternative = "two.sided")
wilcox.test(filter(df_list_dyn_tidy, Sample == "0hr")$Meth, filter(df_list_dyn_tidy, Sample == "24hr")$Meth, alternative = "two.sided")
wilcox.test(filter(df_list_dyn_tidy, Sample == "24hr")$Meth, filter(df_list_dyn_tidy, Sample == "72hr")$Meth, alternative = "two.sided")
wilcox.test(filter(df_list_dyn_tidy, Sample == "24hr")$Meth, filter(df_list_dyn_tidy, Sample == "72hr")$Meth, alternative = "two.sided")


# Build ggplot for static regions 

Boxplot_Meth_StaticRegions <- df_list_stat_tidy %>% 
  ggplot(aes(x = Sample, y = Meth)) + 
  geom_boxplot(color = "black", fill = "grey70", width = 0.50, 
               outlier.alpha = 0.01, outlier.size = 0.75, outlier.shape = NA) + 
  ggtitle("CpG Methylation Among Static Regions") +
  ylab("CpG Methylation Fraction") +
  xlab("Treatment Group (hrs of PMA exposure)") +
  theme_bw() +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=14))

# Write to EPS file 
library(Cairo)
setwd("/Users/kellybarnett/Dropbox (VU Basic Sciences)/Hodges Lab/ATAC_Me_Manuscript_Docs/")
ggsave(filename = "CpG_Methylation_Among_Chromatin_Regions_Openat24hrs_Cluster345_Boxplot.eps",
       plot = print(Boxplot_Meth_OpeningRegions),
       device = cairo_ps)

library(Cairo)
ggsave(filename = "Dropbox (VU Basic Sciences)/Hodges Lab/ATAC_Me_Manuscript_Docs/Mol_Cell_Revisions/CpG_Methylation_Among_Static_Regions_Boxplot.eps",
       plot = print(Boxplot_Meth_StaticRegions),
       device = cairo_ps)


CpG_Methylation_Among_Chromatin_Regions_Openat24hrs_Cluster345_Boxplot.eps
