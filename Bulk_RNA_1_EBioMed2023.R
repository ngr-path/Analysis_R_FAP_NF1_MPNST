#------------------------------------------------------------------------------
# Bulk RNA-seq Analysis (1)
# Title: FAP Expression as a Marker of Malignancy Enabling In-Vivo Imaging 
#        in NF1-Associated Peripheral Nerve Tumors: A Multimodal and Translational Study
#
# Objective:
# Comparative transcriptomic analysis of FAP, GLUT1 (SLC2A1), and GLUT3 (SLC2A3)
# in Neurofibroma vs. MPNST using dataset GSE241224.
#
# Associated Publication in EBioMedicine (2023)
# DOI: 10.1016/j.ebiom.2023.104829
# PMID: 37837931; PMCID: PMC10585232
# Link: https://www.sciencedirect.com/science/article/pii/S235239642300395X
#
# Dataset Source:
# GEO accession: GSE241224
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241224
#
# Author(s): Nic G. Reitsam, Pathology, Faculty of Medicine, University of Augsburg
# Date: 16th June 2025
# ------------------------------------------------------------------------------

library(GEOquery)
library(Biobase)
library(ggplot2)
library(dplyr)
library(broom)
library(patchwork)


#Load GEO dataset
gse <- getGEO("GSE241224", GSEMatrix = TRUE)
gset <- gse[[1]]

#Check if data is already normalized, which is true
summary(exprs(gset))
boxplot(exprs(gset),
        las = 2,           
        outline = FALSE,    
        main = "Sample distributions",
        ylab = "log2(expression)")


feature_data <- fData(gset)
exprs_data <- exprs(gset)
pheno_data <- pData(gset)

#Define relevant probe IDs
#for FAP also see: https://www.thermofisher.com/taqman-gene-expression/product/Hs00990809_m1?CID=&ICID=&subtype=
#for GLUT1 also see: https://www.thermofisher.com/taqman-gene-expression/product/Hs00892681_m1?CID=&ICID=&subtype=
#for GLUT3 also see: https://www.thermofisher.com/order/genome-database/details/gene-expression/Hs00359840_m1
fap_id   <- "TC02002480.hg.1"
glut1_id <- "TC01002578.hg.1"
glut3_id <- "TC12001170.hg.1"

#Expression data
get_expr <- function(probe_id, gene_name) {
  rows <- rownames(exprs_data) %in% probe_id
  vals <- exprs_data[rows, , drop=TRUE]
  tibble(
    gene      = gene_name,
    spot      = colnames(exprs_data),
    expression= as.numeric(vals)
  )
}

expr_df <- bind_rows(
  get_expr(fap_id,   "FAP"),
  get_expr(glut1_id, "SLC2A1"),
  get_expr(glut3_id, "SLC2A3")
)

#condition (MPNST vs Neurofibroma) from pheno_data
expr_df <- expr_df %>%
  mutate(
    condition = ifelse(
      grepl("Neurofibroma", pheno_data$characteristics_ch1.2[ match(spot, rownames(pheno_data)) ]),
      "Neurofibroma", "MPNST"
    ),
    condition = factor(condition, levels = c("Neurofibroma","MPNST")),
    gene      = factor(gene,      levels = c("FAP","SLC2A1","SLC2A3"))
  )

#Plotting function
plot_violin_nm <- function(gene_name) {
  d <- filter(expr_df, gene == gene_name)
  ttest <- t.test(expression ~ condition, data = d)
  pval  <- signif(ttest$p.value, 2)
  
  ggplot(d, aes(x = condition, y = expression, fill = condition)) +
    geom_violin(trim = FALSE, alpha = 0.7, size = 0.3) +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8, color = "black", size = 0.3) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 2.5, color = "white") +
    annotate("text", x = 1.5, y = max(d$expression) * 1.05,
             label = paste0("p=", pval), size = 3, fontface = "bold") +
    scale_fill_manual(values = c("Neurofibroma" = "lightblue","MPNST" = "red")) +
    labs(
      x     = NULL,
      y     = "Expression Level",
      title = gene_name
    ) +
    theme_minimal(base_size = 7) +
    theme(
      plot.title       = element_text(face = "bold", size = 8, hjust = 0.5),
      axis.text.x      = element_text(size = 6),
      axis.text.y      = element_text(size = 6),
      axis.title.y     = element_text(size = 7),
      legend.position  = "none",
      panel.grid.major = element_line(size = 0.2),
      panel.grid.minor = element_blank()
    )
}

# Panels A–C
pA <- plot_violin_nm("FAP")
pB <- plot_violin_nm("SLC2A1") + labs(title = "SLC2A1")
pC <- plot_violin_nm("SLC2A3") + labs(title = "SLC2A3")

#Final Plot
final <- (pA | pB | pC)
theme(plot.tag = element_text(face = "bold", size = 9))

print(final)

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

#gene       n mean_neuro mean_mpnst median_neuro median_mpnst   p_value
#FAP       79       7.73      10.1          7.19        10.6  0.0000234
#SLC2A1    79       9.00       9.44         8.04         8.96 0.432    
#SLC2A3    79      10.1       11.6         10.4         11.9  0.000795 

#FAP expression is markedly higher in MPNST than in benign neurofibromas.
#GLUT1 (SLC2A1) show only a minor, non-significant increase (p=0.432).
#GLUT3 (SLC2A3) is modestly but significantly elevated in MPNSTs - though smaller effect size than FAP.
#FAP-targeted PET may offer superior sensitivity and tumor-to-background contrast in MPNST versus FDG-PET.
