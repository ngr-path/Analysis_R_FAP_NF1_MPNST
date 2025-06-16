#------------------------------------------------------------------------------
# Bulk RNA-seq Analysis (2)
# Title: FAP Expression as a Marker of Malignancy Enabling In-Vivo Imaging 
#        in NF1-Associated Peripheral Nerve Tumors: A Multimodal and Translational Study
#
# Objective:
# Comparative transcriptomic analysis of FAP, GLUT1 (SLC2A1), and GLUT3 (SLC2A3)
# in Neurofibroma vs. MPNST using dataset GSE66743.
#
# Associated Publication in Molecular Oncology 2015
# DOI: 10.1016/j.molonc.2015.02.005
# PMID: 25769404 
# Link: https://febs.onlinelibrary.wiley.com/doi/10.1016/j.molonc.2015.02.005
#
# Dataset Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66743
# GEO accession: GSE66743
# 
#
# Author(s): Nic G. Reitsam, Pathology, Faculty of Medicine, University of Augsburg
# Date: 16th June 2025
# ----------------------------------------------------------------------------------

library(GEOquery)
library(Biobase)
library(ggplot2)
library(dplyr)
library(tibble)
library(patchwork)

# Load data
gse <- getGEO("GSE66743", GSEMatrix = TRUE)
gset <- gse[[1]]

#Check if data is already normalized, which is true
summary(exprs(gset))
boxplot(exprs(gset),
        las = 2,           
        outline = FALSE,    
        main = "Sample distributions",
        ylab = "log2(expression)")

# Expression, features, phenotypes
exprs_data   <- exprs(gset)
feature_data <- fData(gset)
pheno_data   <- pData(gset)

summary(exprs_data)

#Clean condition labels
pheno_data$condition <- ifelse(
  grepl("neurofibroma", tolower(pheno_data$source_name_ch1)), 
  "Neurofibroma", "MPNST"
)
pheno_data$condition <- factor(pheno_data$condition, levels = c("Neurofibroma", "MPNST"))

# Match gene probes
fap_id   <- "hCG1685968.2"
glut1_id <- "hCG40024.3"
glut3_ids <- c("hCG22218.2", "hCG1730826.1")  # use both for GLUT3

# General expression extraction function 
get_expr <- function(gene_ids, gene_name) {
  rows <- feature_data$GeneID %in% gene_ids
  expr_vals <- exprs_data[rows, , drop = FALSE]
  
  if (nrow(expr_vals) == 0) {
    warning(paste("No matching probes found for", gene_name))
    return(tibble(gene = character(), sample = character(), expression = numeric()))
  }
  
  expr_combined <- if (nrow(expr_vals) > 1) {
    colMeans(expr_vals)
  } else {
    as.numeric(expr_vals)
  }
  
  tibble(
    gene = gene_name,
    sample = colnames(exprs_data),
    expression = expr_combined
  )
}

# Special handler for GLUT3 with partial matching
get_expr_glut3 <- function(gene_ids, gene_name) {
  pattern <- paste(gene_ids, collapse = "|")  # regex pattern: "id1|id2"
  rows <- grepl(pattern, feature_data$GeneID)
  expr_vals <- exprs_data[rows, , drop = FALSE]
  
  if (nrow(expr_vals) == 0) {
    warning(paste("No matching probes found for", gene_name))
    return(tibble(gene = character(), sample = character(), expression = numeric()))
  }
  
  expr_combined <- if (nrow(expr_vals) > 1) {
    colMeans(expr_vals)
  } else {
    as.numeric(expr_vals)
  }
  
  tibble(
    gene = gene_name,
    sample = colnames(exprs_data),
    expression = expr_combined
  )
}

# Build combined expression data
expr_df <- bind_rows(
  get_expr(fap_id, "FAP"),
  get_expr(glut1_id, "SLC2A1"),
  get_expr_glut3(glut3_ids, "SLC2A3")
)

# Add condition info
expr_df <- expr_df %>%
  mutate(
    condition = pheno_data$condition[match(sample, rownames(pheno_data))],
    condition = factor(condition, levels = c("Neurofibroma", "MPNST")),
    gene = factor(gene, levels = c("FAP", "SLC2A1", "SLC2A3"))
  )

# Violin plot function
plot_violin_gene <- function(gene_name) {
  d <- filter(expr_df, gene == gene_name)
  t_res <- t.test(expression ~ condition, data = d)
  pval <- signif(t_res$p.value, 2)
  
  ggplot(d, aes(x = condition, y = expression, fill = condition)) +
    geom_violin(trim = FALSE, alpha = 0.7, size = 0.3) +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8, color = "black", size = 0.3) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 2.5, color = "white") +
    annotate("text", x = 1.5, y = max(d$expression) * 1.05,
             label = paste0("p=", pval), size = 3, fontface = "bold") +
    scale_fill_manual(values = c("Neurofibroma" = "lightblue", "MPNST" = "red")) +
    labs(x = NULL, y = "Expression Level", title = gene_name) +
    theme_minimal(base_size = 7) +
    theme(
      plot.title = element_text(face = "bold", size = 8, hjust = 0.5),
      axis.text.x = element_text(size = 6),
      axis.text.y = element_text(size = 6),
      axis.title.y = element_text(size = 7),
      legend.position = "none",
      panel.grid.major = element_line(size = 0.2),
      panel.grid.minor = element_blank()
    )
}

# Generate plots
pA <- plot_violin_gene("FAP")
pB <- plot_violin_gene("SLC2A1")
pC <- plot_violin_gene("SLC2A3")

final_plot <- (pA | pB | pC)
print(final_plot)

#Summary Statistic of Bulk Expression
summary_stats <- expr_df %>%
  group_by(gene) %>%
  group_modify(~ {
    neuro_expr <- .x$expression[.x$condition == "Neurofibroma"]
    mpnst_expr <- .x$expression[.x$condition == "MPNST"]
    t_result   <- t.test(.x$expression ~ .x$condition)
    
    tibble(
      n = nrow(.x),
      mean_neuro = mean(neuro_expr),
      mean_mpnst = mean(mpnst_expr),
      median_neuro = median(neuro_expr),
      median_mpnst = median(mpnst_expr),
      p_value = t_result$p.value
    )
  }) %>%
  ungroup()

print(summary_stats)

# gene       n   mean_neuro  mean_mpnst  median_neuro  median_mpnst   p_value
# FAP        38     -2.25       -0.188       -2.56         0.307      0.00275
# SLC2A1     38     -0.351      -0.112       -0.248        0.200      0.457  
# SLC2A3     38      0.137      -0.0757      -0.0325       0.188      0.726  

# FAP expression is significantly higher in MPNSTs than in neurofibromas (p = 0.0027), suggesting its utility as a malignancy marker.
# SLC2A1 (GLUT1) shows a small, non-significant increase in MPNSTs (p = 0.457).
# SLC2A3 (GLUT3) does not show a meaningful or significant difference between groups (p = 0.726).
# These findings support FAP as a more reliable biomarker and potential imaging target than glucose transporters (FDG-PET/CT) in NF1-associated tumor progression.

