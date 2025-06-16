#Project: FAP as marker for malignancy in NF1 (neurofibroma vs MPNST)
#Spatial transcriptomics data of Juliane Bremer, Pamela Franco, Joelle Aline Menstell, Shelisa Tey, Kamil Kajetan Zajt, Klimentina Popzhelyazkova, Kay Nolte, Jürgen Schlegel, Maria Teresa Pedro, Anja Osterloh, Daniel Delev, Marc Hohenhaus, Christoph Scholz, Oliver Schnell, Juergen Beck, Joachim Weis, Dieter Henrik Heiland, 
#Spatially resolved transcriptomics of benign and malignant peripheral nerve sheath tumors
#Neuro-Oncology, 2025;, noaf016, https://doi.org/10.1093/neuonc/noaf016
#Diagnosis in Supplementary Table S1
#sample1: HPNST, sample2: MPNST, sample3: MPNST, sample 4: HPNST, sample5_ MPNST, sample6: MPNST, sample7: plexiform neurofibroma, sample8: NF2-associated neuropathy

#For First-Pass Translational Validation of FAP expression in MPNST

#In this study, we used spatial transcriptomics data from Bremer et al. to evaluate the expression of Fibroblast Activation Protein (FAP) in malignant peripheral nerve sheath tumors (MPNST) and its relationship to tumor cell markers, including glucose transporters relevant for FDG-PET imaging (e.g., SLC2A1, SLC2A3). 
#Rather than applying advanced spatial modeling techniques, we employed a deliberately straightforward and interpretable strategy.
#Co-expression analysis between FAP and MPNST-associated genes was performed using Fisher’s exact test and Spearman correlation, focusing on the biological relevance of shared expression.
#Joint expression fractions (spots positive for both FAP and another marker) were computed to approximate spatial overlap and potential co-localization.
#This minimal statistical framework ensures robust and reproducible interpretation while highlighting biologically meaningful expression patterns of FAP in relation to tumor cells and glucose metabolism.
#FAP and other marker expression
#what to report: Marker /	Nonzero Spots /	Fisher P-value /	Spearman ρ / Joint+ (%)
# Author(s): Nic G. Reitsam, Pathology, Faculty of Medicine, University of Augsburg
# Date: 16th June 2025

#installing packages
install.packages("BiocManager")
BiocManager::install("Seurat")
BiocManager::install("Matrix")
install.packages("ggplot2")
install.packages("patchwork") 

library(patchwork)
library(Seurat)
library(ggplot2)
library(Matrix)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggplot2)
library(viridis)
library(patchwork)
library(grid)


# Sample 1 (HPNST)
data_dir1 <- "/Users/nicreitsam/Downloads/SpatialTranscriptomics/sample1/"
sample1 <- Load10X_Spatial(data.dir = data_dir1)
sample1 <- subset(sample1, subset = nFeature_Spatial > 50 & nCount_Spatial < 20000)
sample1 <- NormalizeData(sample1, normalization.method = "LogNormalize", scale.factor = 10000)
sample1 <- FindVariableFeatures(sample1, selection.method = "vst", nfeatures = 2000)
sample1 <- ScaleData(sample1)
SpatialFeaturePlot(sample1, features = "FAP", pt.size.factor = 1.6) + ggtitle("Sample1 (HPNST): FAP")

# Sample 2 (MPNST)
data_dir2 <- "/Users/nicreitsam/Downloads/SpatialTranscriptomics/sample2/"
sample2 <- Load10X_Spatial(data.dir = data_dir2)
sample2 <- subset(sample2, subset = nFeature_Spatial > 50 & nCount_Spatial < 20000)
sample2 <- NormalizeData(sample2, normalization.method = "LogNormalize", scale.factor = 10000)
sample2 <- FindVariableFeatures(sample2, selection.method = "vst", nfeatures = 2000)
sample2 <- ScaleData(sample2)
SpatialFeaturePlot(sample2, features = "FAP", pt.size.factor = 1.6) + ggtitle("Sample2 (MPNST): FAP")

# Sample 3 (MPNST)
data_dir3 <- "/Users/nicreitsam/Downloads/SpatialTranscriptomics/sample3/"
sample3 <- Load10X_Spatial(data.dir = data_dir3)
sample3 <- subset(sample3, subset = nFeature_Spatial > 50 & nCount_Spatial < 20000)
sample3 <- NormalizeData(sample3, normalization.method = "LogNormalize", scale.factor = 10000)
sample3 <- FindVariableFeatures(sample3, selection.method = "vst", nfeatures = 2000)
sample3 <- ScaleData(sample3)
SpatialFeaturePlot(sample3, features = "FAP", pt.size.factor = 1.6) + ggtitle("Sample3 (MPNST): FAP")

# Sample 4 (HPNST)
data_dir4 <- "/Users/nicreitsam/Downloads/SpatialTranscriptomics/sample4/"
sample4 <- Load10X_Spatial(data.dir = data_dir4)
sample4 <- subset(sample4, subset = nFeature_Spatial > 50 & nCount_Spatial < 20000)
sample4 <- NormalizeData(sample4, normalization.method = "LogNormalize", scale.factor = 10000)
sample4 <- FindVariableFeatures(sample4, selection.method = "vst", nfeatures = 2000)
sample4 <- ScaleData(sample4)
SpatialFeaturePlot(sample4, features = "FAP", pt.size.factor = 1.6) + ggtitle("Sample4 (HPNST): FAP")

# Sample 5 (MPNST)
data_dir5 <- "/Users/nicreitsam/Downloads/SpatialTranscriptomics/sample5/"
sample5 <- Load10X_Spatial(data.dir = data_dir5)
sample5 <- subset(sample5, subset = nFeature_Spatial > 50 & nCount_Spatial < 20000)
sample5 <- NormalizeData(sample5, normalization.method = "LogNormalize", scale.factor = 10000)
sample5 <- FindVariableFeatures(sample5, selection.method = "vst", nfeatures = 2000)
sample5 <- ScaleData(sample5)
SpatialFeaturePlot(sample5, features = "FAP", pt.size.factor = 1.6) + ggtitle("Sample5 (MPNST): FAP")

# Sample 6 (MPNST)
data_dir6 <- "/Users/nicreitsam/Downloads/SpatialTranscriptomics/sample6/"
sample6 <- Load10X_Spatial(data.dir = data_dir6)
sample6 <- subset(sample6, subset = nFeature_Spatial > 50 & nCount_Spatial < 20000)
sample6 <- NormalizeData(sample6, normalization.method = "LogNormalize", scale.factor = 10000)
sample6 <- FindVariableFeatures(sample6, selection.method = "vst", nfeatures = 2000)
sample6 <- ScaleData(sample6)
SpatialFeaturePlot(sample6, features = "FAP", pt.size.factor = 1.6) + ggtitle("Sample6 (MPNST): FAP")

# Sample 7 (plexiform neurofibroma)
data_dir7 <- "/Users/nicreitsam/Downloads/SpatialTranscriptomics/sample7/"
sample7 <- Load10X_Spatial(data.dir = data_dir7)
sample7 <- subset(sample7, subset = nFeature_Spatial > 50 & nCount_Spatial < 20000)
sample7 <- NormalizeData(sample7, normalization.method = "LogNormalize", scale.factor = 10000)
sample7 <- FindVariableFeatures(sample7, selection.method = "vst", nfeatures = 2000)
sample7 <- ScaleData(sample7)
SpatialFeaturePlot(sample7, features = "FAP", pt.size.factor = 1.6) + ggtitle("Sample7 (Plexiform NF): FAP")

# Sample 8 (NF2-associated neuropathy)
data_dir8 <- "/Users/nicreitsam/Downloads/SpatialTranscriptomics/sample8/"
sample8 <- Load10X_Spatial(data.dir = data_dir8)
sample8 <- subset(sample8, subset = nFeature_Spatial > 50 & nCount_Spatial < 20000)
sample8 <- NormalizeData(sample8, normalization.method = "LogNormalize", scale.factor = 10000)
sample8 <- FindVariableFeatures(sample8, selection.method = "vst", nfeatures = 2000)
sample8 <- ScaleData(sample8)
SpatialFeaturePlot(sample8, features = "FAP", pt.size.factor = 1.6) + ggtitle("Sample8 (NF2 Neuropathy): FAP")


sample1$lesion_type <- "HPNST"
sample2$lesion_type <- "MPNST"
sample3$lesion_type <- "MPNST"
sample4$lesion_type <- "HPNST"
sample5$lesion_type <- "MPNST"
sample6$lesion_type <- "MPNST"
sample7$lesion_type <- "plexiform_neurofibroma"
sample8$lesion_type <- "NF2_associated_neuropathy"

# Merge all samples into  Seurat object
combined <- merge(
  sample1, y = list(sample2, sample3, sample4, sample5, sample6, sample7, sample8),
  add.cell.ids = paste0("s", 1:8), project = "Visium_NF1"
)

#"high‑FAP" as those spots with FAP expression above the 75th percentile
fap_values <- FetchData(combined, vars = "FAP")$FAP
threshold <- quantile(fap_values, 0.75)
combined$FAP_high <- fap_values > threshold

#Contingency table & Fisher’s exact test
contab <- table(
  FAP_high = combined$FAP_high,
  Malignant = combined$lesion_type == "MPNST"
)
print(contab)

fisher_res <- fisher.test(contab)
print(fisher_res)

prop_df <- combined@meta.data %>%
  as_tibble(rownames = "spot") %>%
  mutate(Malignant = lesion_type == "MPNST") %>%
  group_by(Malignant) %>%
  summarize(
    total_spots   = n(),
    highFAP_spots = sum(FAP_high),
    prop_highFAP  = highFAP_spots / total_spots * 100
  ) %>%
  mutate(
    Malignant = if_else(Malignant, "MPNST", "Non-MPNST")
  )

#Bar plot of % high-FAP spots
ggplot(prop_df, aes(x = Malignant, y = prop_highFAP, fill = Malignant)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = paste0(round(prop_highFAP,1), "%")),
            vjust = -0.5, size = 5) +
  scale_y_continuous(expand = expansion(c(0,0.1))) +
  labs(
    x = NULL,
    y = "Percentage of High-FAP Spots",
    title = "High-FAP Spot Frequency\nMPNST vs Non-MPNST"
  ) +
  theme_minimal(base_size = 14)

#Annotate with odds ratio & p-value
or <- round(fisher_res$estimate, 2)
pv <- signif(fisher_res$p.value, 2)
ggplot(prop_df, aes(x = Malignant, y = prop_highFAP, fill = Malignant)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = paste0(round(prop_highFAP,1), "%")),
            vjust = -0.5, size = 5) +
  annotate(
    "text", x = 2, y = max(prop_df$prop_highFAP) + 5,
    label = paste0("OR=", or, ", p=", pv),
    hjust = 1, size = 5
  ) +
  scale_y_continuous(expand = expansion(c(0,0.2))) +
  labs(
    x = NULL,
    y = "Percentage of High-FAP Spots",
    title = "High-FAP Spot Frequency\nMPNST vs Non-MPNST"
  ) +
  theme_minimal(base_size = 14)

# Total spots
total_spots <- ncol(combined)
cat("Total spots:", total_spots, "\n")

#Spots by lesion type
spots_by_lesion <- table(combined$lesion_type)
print(spots_by_lesion)

#malignant vs. non‐malignant (MPNST vs others):
combined$malignant <- combined$lesion_type == "MPNST"
spots_by_malignancy <- table(combined$malignant)
names(spots_by_malignancy) <- c("Non-MPNST", "MPNST")
print(spots_by_malignancy)

#Sum of non-MPNST and MPNST spots
cat("Non-MPNST spots:", spots_by_malignancy["Non-MPNST"], "\n")
cat("MPNST spots:",     spots_by_malignancy["MPNST"],     "\n")

genes_to_test <- c("FAP", "SLC2A1", "SLC2A3")
results <- list()

for (gene in genes_to_test) {
  expr <- FetchData(combined, vars = gene)[[gene]]
  threshold <- quantile(expr, 0.75, na.rm = TRUE)
  high_expr <- expr > threshold
  
  # Save into meta.data
  combined[[paste0(gene, "_high")]] <- high_expr
  
  # Contingency table and Fisher's test
  contab <- table(
    high = high_expr,
    Malignant = combined$lesion_type == "MPNST"
  )
  
  fisher_res <- fisher.test(contab)
  
  results[[gene]] <- list(
    table = contab,
    OR = fisher_res$estimate,
    p = fisher_res$p.value
  )
}

# how summary
for (gene in genes_to_test) {
  cat(paste0(
    "\n--- ", gene, " ---\n",
    "Odds Ratio: ", round(results[[gene]]$OR, 2), "\n",
    "p-value: ", signif(results[[gene]]$p, 3), "\n"
  ))
}


#Generalized plotting function with y-axis limited to 20%
plot_marker_enrichment <- function(marker, combined_obj, fisher_results) {
  high_col <- paste0(marker, "_high")
  
  prop_df <- combined_obj@meta.data %>%
    as_tibble(rownames = "spot") %>%
    mutate(Malignant = lesion_type == "MPNST") %>%
    group_by(Malignant) %>%
    summarize(
      total_spots   = n(),
      high_spots    = sum(.data[[high_col]]),
      prop_high     = high_spots / total_spots * 100
    ) %>%
    mutate(
      Malignant = if_else(Malignant, "MPNST", "Non-MPNST")
    )
  
  or <- round(fisher_results[[marker]]$OR, 2)
  pv <- signif(fisher_results[[marker]]$p, 2)
  
  ggplot(prop_df, aes(x = Malignant, y = prop_high, fill = Malignant)) +
    geom_col(width = 0.6, show.legend = FALSE) +
    geom_text(aes(label = paste0(round(prop_high, 1), "%")),
              vjust = -0.5, size = 5) +
    annotate(
      "text", x = 2, y = 19,
      label = paste0("OR=", or, ", p=", pv),
      hjust = 1, size = 5
    ) +
    scale_y_continuous(limits = c(0, 20), expand = c(0, 0)) +
    labs(
      x = NULL,
      y = paste0("Percentage of High-", marker, " Spots"),
      title = paste0("High-", marker, " Spot Frequency\nMPNST vs Non-MPNST")
    ) +
    theme_minimal(base_size = 14)
}


plot_marker_enrichment("FAP", combined, results)
plot_marker_enrichment("SLC2A1", combined, results)
plot_marker_enrichment("SLC2A3", combined, results)

# List of MPNST samples
mpnst_samples <- list(
  sample2 = sample2,
  sample3 = sample3,
  sample5 = sample5,
  sample6 = sample6
)

# Markers for tumor states and stroma
markers <- c("APOD", "NGFR", "L1CAM", "S100B", "VIM", "MYH9", "ZEB1", "CDK4", "DCN", "LAMA2", "COL1A1", "PDGFRB", "SOX10")

#Initialize empty result table
all_results <- data.frame()

#Loop through MPNST samples
for (sample_name in names(mpnst_samples)) {
  sample <- mpnst_samples[[sample_name]]
  

  expr_df <- FetchData(sample, vars = c("FAP", markers))
  

  for (gene in markers) {
    fap_bin <- expr_df$FAP > 0
    gene_bin <- expr_df[[gene]] > 0
    
    if (sum(gene_bin) > 0) {
      fisher_p <- fisher.test(table(fap_bin, gene_bin))$p.value
      rho <- cor(expr_df$FAP, expr_df[[gene]], method = "spearman")
      joint_pos <- sum(fap_bin & gene_bin) / sum(fap_bin) * 100
      nonzero_spots <- sum(gene_bin)
      
      all_results <- rbind(all_results, data.frame(
        Sample = sample_name,
        Marker = gene,
        Nonzero_Spots = nonzero_spots,
        Fisher_P = signif(fisher_p, 3),
        Spearman_rho = round(rho, 2),
        Joint_Positive_Perc = round(joint_pos, 1)
      ))
    }
  }
}

print(all_results)


all_results %>%
  mutate(log_p = -log10(as.numeric(Fisher_P))) %>%
  ggplot(aes(x = Marker, y = Sample, size = Joint_Positive_Perc, color = log_p)) +
  geom_point() +
  scale_color_gradient(low = "lightgrey", high = "darkred") +
  scale_size(range = c(3,10)) +
  theme_minimal(base_size = 13) +
  labs(title = "FAP Co-Expression: Joint Positivity & Significance",
       x = "Marker", y = "MPNST Sample", color = "-log10(p)", size = "% Joint+ Spots") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Add a column for marker category (Tumor vs Stroma)
all_results <- all_results %>%
  mutate(Category = case_when(
    Marker %in% c("APOD", "NGFR", "L1CAM", "S100B", "VIM", "MYH9", "ZEB1", "CDK4", "SOX10") ~ "Tumor Cell Marker",
    Marker %in% c("DCN", "LAMA2", "COL1A1", "PDGFRB") ~ "Stroma Marker",
    TRUE ~ "Other"
  ))

all_results %>%
  mutate(log_p = -log10(as.numeric(Fisher_P))) %>%
  ggplot(aes(x = Marker, y = Sample, size = Joint_Positive_Perc, color = log_p)) +
  geom_point() +
  scale_color_gradient(low = "lightblue", high = "red") +
  scale_size(range = c(3, 10)) +
  theme_minimal(base_size = 13) +
  facet_wrap(~Category, scales = "free_x") +
  labs(
    title = "FAP Co-Expression: Tumor vs Stroma Markers",
    x = "Marker", y = "MPNST Sample",
    color = "-log10(p)", size = "% Joint+ Spots"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Palette
cb_palette <- c("Non-MPNST" = "lightblue",   # blue
                "MPNST"     = "red")   # vermillion

#Plotting function
plot_marker_enrichment_nm <- function(marker, combined_obj, fisher_results) {
  high_col <- paste0(marker, "_high")
  
  prop_df <- combined_obj@meta.data %>%
    as_tibble(rownames = "spot") %>%
    mutate(Malignant = if_else(lesion_type == "MPNST", "MPNST", "Non-MPNST")) %>%
    group_by(Malignant) %>%
    summarise(
      n_spots    = n(),
      n_high     = sum(.data[[high_col]]),
      pct_high   = n_high / n_spots * 100
    ) %>%
    ungroup()
  
  or  <- round(fisher_results[[marker]]$OR, 2)
  pv  <- signif(fisher_results[[marker]]$p, 2)
  anno_text <- paste0("OR=", or, "  p=", pv)
  
  ggplot(prop_df, aes(x = Malignant, y = pct_high, fill = Malignant)) +
    geom_col(width = 0.5, show.legend = FALSE) +
    geom_text(aes(label = paste0(round(pct_high,1), "%")), 
              vjust = -0.5, size = 3.5, color = "black") +
    annotate("text", x = 2, y = 19.5, label = anno_text, size = 3.5, hjust = 1) +
    scale_fill_manual(values = cb_palette) +
    scale_y_continuous(limits = c(0, 20), expand = c(0, 0)) +
    labs(
      x = NULL,
      y = "High-Expression Spots (%)",
      title = paste0(marker, " high-expression frequency")
    ) +
    theme_minimal(base_size = 8) +
    theme(
      plot.title       = element_text(face = "bold", size = 9, hjust = 0.5),
      axis.text.x      = element_text(size = 8),
      axis.title.y     = element_text(size = 8),
      panel.grid.major = element_line(size = 0.2),
      panel.grid.minor = element_blank()
    )
}

pA <- plot_marker_enrichment_nm("FAP",   combined, results)
pB <- plot_marker_enrichment_nm("SLC2A1", combined, results)
pC <- plot_marker_enrichment_nm("SLC2A3", combined, results)

final_figure <- (pA | pB | pC)

grid.newpage()
grid.draw(final_figure)


p_coexpr_filled <- all_results %>%
  mutate(log_p = -log10(as.numeric(Fisher_P))) %>%
  ggplot(aes(x = Marker, y = Sample, size = Joint_Positive_Perc)) +
  
  geom_point(aes(fill = log_p), shape = 21, color = "black", stroke = 0.2) +
  scale_fill_viridis(option = "D", direction = -1, name = expression(-log[10](p))) +
  scale_size_continuous(range = c(2, 8), name = "% Joint-positive") +
  facet_wrap(~ Category, scales = "free_x", strip.position = "bottom") +
  labs(
    x = "Marker",
    y = "MPNST samples",
    title = "FAP co-expression with tumor vs stromal markers"
  ) +
  theme_minimal(base_size = 7) +
  theme(
    plot.title         = element_text(face = "bold", size = 8, hjust = 0.5),
    axis.title         = element_text(size = 8),
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y        = element_text(size = 6),
    strip.placement    = "outside",
    strip.text         = element_text(size = 7, face = "bold"),
    legend.title       = element_text(size = 7),
    legend.text        = element_text(size = 6),
    panel.grid.major   = element_line(size = 0.2),
    panel.grid.minor   = element_blank()
  )

p_coexpr_filled + 
  theme(plot.tag = element_text(face = "bold", size = 8))


#Methods:
#We analyzed spatial transcriptomics (10x Visium) data from eight NF1-associated lesions (two HPNST, four MPNST, one plexiform neurofibroma, one NF2-associated neuropathy) using Seurat (v4). 
#After filtering to spots with >50 genes and <20 000 UMIs, we log-normalized (scale.factor=10 000), identified 2 000 variable features (vst), and scaled all genes. 
#We then dichotomized FAP expression at its 75th percentile to define “high-FAP” spots and constructed a 2×2 contingency table of high-FAP versus MPNST status. 
#Enrichment of FAP-high spots in MPNST regions was assessed by Fisher’s exact test (odds ratio, 95% CI, p-value). 
#To visualize this, we plotted the percentage of high-FAP spots in malignant versus non-malignant lesions and overlaid high-FAP spots on the spatial tissue map. 
#This minimal, reproducible framework highlights FAP’s potential as a specific marker of malignant transformation in NF1.


