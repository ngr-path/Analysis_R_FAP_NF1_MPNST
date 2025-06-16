#TCGA Bulk RNA Data of Different Sarcoma Types (5)
#Title: FAP Expression as a Marker of Malignancy Enabling In-Vivo Imaging 
#        in NF1-Associated Peripheral Nerve Tumors: A Multimodal and Translational Study
#
# Publication: doi: 10.1016/j.cell.2017.10.014 in Cell 2017: TCGA Sarcoma
# Data from: https://www.cbioportal.org/study/summary?id=sarc_tcga_pan_can_atlas_2018
# Batch-normalized RSEM RNA-seq data 
# - Data wrangling and visualization of FAP mRNA expression
# - Statistical analysis using Kruskal-Wallis and Dunn's test
# Author(s): Nic G. Reitsam, Pathology, Faculty of Medicine, University of Augsburg
# Date: 16th June 2025

library(dplyr)
library(ggplot2)
library(ggpubr)
library(FSA)

# Load and merge data
cancer_df <- read.delim("/Users/nicreitsam/Downloads/Cancer_Type_Detailed.txt", header = TRUE, stringsAsFactors = FALSE)
expr_df <- read.delim("/Users/nicreitsam/Downloads/FAP__mRNA_Expression,_RSEM_(Batch_normalized_from_Illumina_HiSeq_RNASeqV2).txt", header = TRUE, stringsAsFactors = FALSE)

# expression data with cancer type
merged_df <- merge(expr_df, cancer_df, by = "Sample.ID")
colnames(merged_df)[which(grepl("FAP", colnames(merged_df)))] <- "FAP_expression"

# shortened sarcoma subtype labels
merged_df <- merged_df %>%
  mutate(Sarcoma_Short = recode(Cancer.Type.Detailed,
                                "Dedifferentiated Liposarcoma" = "DDLPS",
                                "Leiomyosarcoma" = "LMS",
                                "Myxofibrosarcoma" = "MFS",
                                "Undifferentiated Pleomorphic Sarcoma/Malignant Fibrous Histiocytoma/High-Grade Spindle Cell Sarcoma" = "UPS/MFH",
                                "Synovial Sarcoma" = "SS",
                                "Malignant Peripheral Nerve Sheath Tumor" = "MPNST",
                                "Desmoid/Aggressive Fibromatosis" = "Desmoid",
                                .default = "Other"
  ))

# Log2-transform FAP expression
merged_df$log_FAP_expression <- log2(merged_df$FAP_expression + 1)

# Calculate sample counts
sample_counts <- merged_df %>%
  group_by(Sarcoma_Short) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = paste0(Sarcoma_Short, "\n(n=", n, ")"))

# Merge sample counts back into the main dataset
merged_df <- merged_df %>%
  left_join(sample_counts, by = "Sarcoma_Short")

# Create the boxplot/jitterplot
ggplot(merged_df, aes(x = label, y = log_FAP_expression)) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.5, color = "gray30") +
  geom_boxplot(width = 0.2, outlier.shape = NA, fill = "white", color = "black") +
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 3, color = "black") +
  labs(
    title = "FAP Expression Across Sarcoma Subtypes in TCGA",
    x = NULL,
    y = expression(log[2]*"(FAP Expression + 1)")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.title.y = element_text(size = 13),
    plot.title = element_text(size = 16, face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  stat_compare_means(method = "kruskal.test", label.y = max(merged_df$log_FAP_expression, na.rm = TRUE) + 1.5,
                     label = "p.format", size = 5)

# Perform Kruskal-Wallis test for overall significance
kruskal_test_result <- kruskal.test(log_FAP_expression ~ Sarcoma_Short, data = merged_df)
print(paste("Kruskal-Wallis test p-value:", format(kruskal_test_result$p.value, digits = 3)))

# Perform pairwise comparisons using Dunn's test (Bonferroni adjustment)
pairwise_comparisons <- dunnTest(log_FAP_expression ~ Sarcoma_Short, data = merged_df, method = "bonferroni")
print(pairwise_comparisons)

#The Kruskal-Wallis test revealed a highly significant overall difference in FAP expression across sarcoma subtypes.
#MPNSTs show relevant FAP mRNA expression

sessionInfo()