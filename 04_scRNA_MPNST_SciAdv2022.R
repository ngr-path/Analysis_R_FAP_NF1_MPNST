#Single-Cell Data Analysis (4)
#Title: FAP Expression as a Marker of Malignancy Enabling In-Vivo Imaging 
#        in NF1-Associated Peripheral Nerve Tumors: A Multimodal and Translational Study
#
#Associated publication: Wu LMN, Zhang F, Rao R, Adam M et al. Single-cell multiomics identifies clinically relevant mesenchymal stem-like cells and key regulators for MPNST malignancy. Sci Adv 2022 Nov 4;8(44):eabo5442. PMID: 36322658
#To investigate the cellular transition from human benign NF to MPNST, we generated single-cell transcriptomes from human NF1-mutated PNF (n = 10; 55,770 cells)
#and NF1-associated MPNST (n = 4; 22,661 cells) samples (Fig. 5A, fig. S7A, and data S3)

#Aim: Studying FAP expression on NF1-associated MPNST vs neurofibroma on single cell level
#
#Dataset Source:
#GEO accession: GSE178989
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178989
#
#Further Information
#Suppiah et al, Nature Communications, Two MPNST subtypes: https://doi.org/10.1038/s41467-023-38432-6
#MPNST-G1: early neural crest cell specification (TWIST1, SOX9, SNAI2, OTX2, PAX3 and PAX6)
#MPNST-G2: Schwann-cell precursor-like cell phenotype by overexpressing GAP43, PLP1 and NGFR, with an absence of Schwann cell markers (S100B, PMP2, ERBB3, MPZ)

# Author(s): Nic G. Reitsam, Pathology, Faculty of Medicine, University of Augsburg
# Date: 16th June 2025

library(Seurat)
library(glmGamPoi)
library(purrr)
library(patchwork)
library(ggplot2)
library(slingshot)

base_dir <- "/Users/nicreitsam/Downloads/GSE179033"

sample_folders <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

#  14 samples (GSM5404031 to GSM5404044)
sample_folders <- sample_folders[grepl("GSM54040(3[1-9]|4[0-4])", basename(sample_folders))]

sample_folders <- sample_folders[order(sample_folders)]

#List Seurat objects
seurat_list <- list()

for (i in seq_along(sample_folders)) {
  folder <- sample_folders[i]
  sample_name <- basename(folder)
  message("Loading sample ", i, ": ", sample_name)
  
  mat <- Read10X(data.dir = folder)
  seurat_obj <- CreateSeuratObject(counts = mat, project = sample_name)
  
  # Assign tumor type based on sample name
  if (grepl("NF", sample_name, ignore.case = TRUE)) {
    seurat_obj$tumor_type <- "NF"
  } else if (grepl("MPNST", sample_name, ignore.case = TRUE)) {
    seurat_obj$tumor_type <- "MPNST"
  } else {
    seurat_obj$tumor_type <- "Unknown"
  }
  
  seurat_obj$sample <- sample_name
  
  seurat_list[[paste0("Seurat", i)]] <- seurat_obj
}

for(i in seq_along(seurat_list)) {
  seurat_list[[i]][["percent.mt"]] <- PercentageFeatureSet(seurat_list[[i]], pattern = "^MT-")
}

for(i in seq_along(seurat_list)) {
  seurat_list[[i]] <- subset(seurat_list[[i]],
                             subset = nFeature_RNA > 200 &
                               nFeature_RNA < 9000 &
                               percent.mt < 20)
}

for (i in seq_along(seurat_list)) {
  seurat_list[[i]] <- NormalizeData(seurat_list[[i]])
  seurat_list[[i]] <- FindVariableFeatures(seurat_list[[i]])
  seurat_list[[i]] <- ScaleData(seurat_list[[i]], vars.to.regress = "percent.mt")
  seurat_list[[i]] <- RunPCA(seurat_list[[i]], verbose = FALSE)
}

#tumor cell markers
tumor_markers <- c("APOD", "NGFR", "L1CAM", "S100B", "VIM", "MYH9", "ZEB1", "CDK4", "SOX10")

#only MPNST samples
mpnst_list <- lapply(seurat_list, function(x) if (x$tumor_type[1] == "MPNST") x else NULL)
mpnst_list <- Filter(Negate(is.null), mpnst_list)

# Merge MPNST samples
mpnst_combined <- mpnst_list[[1]]
if (length(mpnst_list) > 1) {
  for (i in 2:length(mpnst_list)) {
    mpnst_combined <- merge(mpnst_combined, y = mpnst_list[[i]])
    gc()
  }
}


mpnst_combined <- NormalizeData(mpnst_combined, verbose = FALSE)

mpnst_combined <- FindVariableFeatures(mpnst_combined, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

if ("percent.mt" %in% colnames(mpnst_combined@meta.data)) {
  mpnst_combined <- ScaleData(mpnst_combined, vars.to.regress = "percent.mt", verbose = FALSE)
} else {
  mpnst_combined <- ScaleData(mpnst_combined, verbose = FALSE)
}

mpnst_combined <- RunPCA(mpnst_combined, npcs = 30, verbose = FALSE)

mpnst_combined <- RunUMAP(mpnst_combined, dims = 1:20, verbose = FALSE)

FeaturePlot(mpnst_combined, features = "FAP", cols = c("lightgrey", "red"), reduction = "umap") +
  ggtitle("FAP expression in MPNST samples")


mpnst_combined <- FindNeighbors(mpnst_combined, dims = 1:20)
mpnst_combined <- FindClusters(mpnst_combined, resolution = 0.1)  # Adjust resolution for more/fewer clusters
mpnst_combined <- JoinLayers(mpnst_combined)

DimPlot(mpnst_combined, reduction = "umap", label = TRUE) + ggtitle("Clusters in MPNST samples")

mpnst_combined

cluster_markers <- FindAllMarkers(
  mpnst_combined,
  only.pos = TRUE,         # Only positive markers
  min.pct = 0.25,          # Expressed in at least 25% of cells
  logfc.threshold = 0.25   # Log fold-change cutoff
)

top10_markers <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

print(top10_markers, n=130)

cluster_ids <- c(
  "0" = "Tumor_Mesenchymal_like",
  "1" = "Macrophages_TAMs",
  "2" = "Tumor_Schwann_like",
  "3" = "Tumor_Neurodevelopmental",
  "4" = "T_cells",
  "5" = "Fibroblasts",
  "6" = "Endothelial",
  "7" = "Monocytes",
  "8" = "B_cells",
  "9" = "Pericytes_or_MesenchymalStromal",
  "10" = "Cycling_Tumor",
  "11" = "Smooth_muscle_or_Myofibroblasts",
  "12" = "pDCs_or_NK_like"
)


mpnst_combined <- RenameIdents(mpnst_combined, cluster_ids)

levels(Idents(mpnst_combined))

umap_celltype <- DimPlot(mpnst_combined, label = TRUE, repel = TRUE, 
                         label.size = 4, group.by = "ident") + 
  ggtitle("Cell Types") + 
  theme(plot.title = element_text(size = 18))

#FAP expression
umap_fap <- FeaturePlot(mpnst_combined, features = "FAP", 
                        order = TRUE, cols = c("lightgray", "red")) + 
  ggtitle("FAP Expression") + 
  theme(plot.title = element_text(size = 18))

#SLC2A1 expression
umap_slc2a1 <- FeaturePlot(mpnst_combined, features = "SLC2A1", 
                           order = TRUE, cols = c("lightgray", "red")) + 
  ggtitle("SLC2A1 Expression") + 
  theme(plot.title = element_text(size = 18))

#SLC2A3 expression
umap_slc2a3 <- FeaturePlot(mpnst_combined, features = "SLC2A3", 
                           order = TRUE, cols = c("lightgray", "red")) + 
  ggtitle("SLC2A3 Expression") + 
  theme(plot.title = element_text(size = 18))

umap_celltype
umap_fap
umap_slc2a1
umap_slc2a3



Idents(mpnst_combined) <- factor(Idents(mpnst_combined), levels = names(cluster_ids))

# Dot plot for FAP, SLC2A1, SLC2A3
DotPlot(mpnst_combined, features = c("FAP", "SLC2A1", "SLC2A3")) + 
  RotatedAxis() +
  labs(title = "Marker Expression") +
  theme(plot.title = element_text(size = 14))

#Violin plot for FAP
VlnPlot(mpnst_combined, features = "FAP")

pan_tumor <- c("APOD", "NGFR", "S100B", "VIM", "MYH9", "CDK4")
cluster_markers_simple <- list(
  Tumor_Mesenchymal_like        = "ZEB1",           # or ZEB1
  Macrophages_TAMs              = "CD68",
  Schwann_like                  = "SOX10",
  Tumor_Neurodevelopmental      = "L1CAM",         # or APOD
  T_cells                       = "CD3D",
  Fibroblasts                   = "DCN","COL1A1",           # Decorin
  Endothelial                   = "CD34",        # CD31
  Monocytes                     = "CD14","LYZ",
  B_cells                       = "CD79A",
  Pericytes_or_MesenchymalStromal = "PDGFRB",
  Cycling_Tumor                = "MKI67",
  Smooth_muscle_or_Myofibroblasts = "ACTA2",
  pDCs_or_NK_like               = "GNLY", "NKG7", "CLEC4C", "PRF1"
)


# Combine for plotting
marker_map <- c(
  setNames(as.list(pan_tumor), paste0("Tumor_Marker_", pan_tumor)),
  cluster_markers_simple
)

plots <- imap(marker_map, function(genes, name) {
  FeaturePlot(
    mpnst_combined,
    features = genes,
    cols = c("lightgrey", "darkred"),
    reduction = "umap",
    combine = FALSE
  ) %>%
    wrap_plots(ncol = length(genes)) &
    plot_annotation(title = name) &
    theme(plot.title = element_text(size = 14))
})


# Show plots 1–4
wrap_plots(plots[1:4], ncol = 2)

# Show plots 5–8
wrap_plots(plots[5:8], ncol = 2)

# Show plots 9–12
wrap_plots(plots[9:12], ncol = 2)

# Show remaining
wrap_plots(plots[13:16], ncol = 2)
wrap_plots(plots[17:20], ncol = 2)
wrap_plots(plots[21:24], ncol = 2)

tumor_clusters <- c("Tumor_Mesenchymal_like", "Tumor_Neurodevelopmental", "Cycling_Tumor")

avg_expr <- AverageExpression(mpnst_combined, features = c("FAP", "SLC2A1", "SLC2A3"))$RNA
avg_expr <- as.data.frame(t(avg_expr))
avg_expr$Marker <- rownames(avg_expr)
avg_expr_long <- tidyr::pivot_longer(avg_expr, cols = -Marker, names_to = "Cluster", values_to = "Expression")

ggplot(avg_expr_long, aes(x = Cluster, y = Expression, fill = Marker)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Average Expression per Cluster", y = "Avg Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Cluster identities
cluster_counts <- table(Idents(mpnst_combined))

# Convert to data.frame for ggplot
df_counts <- as.data.frame(cluster_counts)
colnames(df_counts) <- c("Cluster", "CellCount")

# Barplot
ggplot(df_counts, aes(x = Cluster, y = CellCount, fill = Cluster)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Number of Cells per Cluster in MPNST Data",
       x = "Cluster",
       y = "Number of Cells") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1),
        legend.position = "none")

mpnst_combined

print(df_counts)

#total cells: 22664
#print(df_counts)
#Cluster CellCount
#1           Tumor_Mesenchymal_like      5134
#2                 Macrophages_TAMs      3451
#3                     Schwann_like      3123
#4         Tumor_Neurodevelopmental      2565
#5                          T_cells      2006
#6                      Fibroblasts      1444
#7                      Endothelial      1141
#8                        Monocytes      1017
#9                          B_cells       923
#10 Pericytes_or_MesenchymalStromal       670
#11                   Cycling_Tumor       552
#12 Smooth_muscle_or_Myofibroblasts       482
#13                 pDCs_or_NK_like       156

#Threshold for FAP-positivity
fap_positive_threshold <- 0.5

# Add metadata column for FAP+ cells
mpnst_combined$FAP_positive <- FetchData(mpnst_combined, "FAP") > fap_positive_threshold

# % of FAP+ cells per cluster
fap_prop <- as.data.frame(table(Idents(mpnst_combined), mpnst_combined$FAP_positive))
colnames(fap_prop) <- c("Cluster", "FAP_Positive", "Count")
fap_prop <- tidyr::pivot_wider(fap_prop, names_from = FAP_Positive, values_from = Count, values_fill = 0)
colnames(fap_prop) <- c("Cluster", "FAP_Negative", "FAP_Positive")
fap_prop$Total <- fap_prop$FAP_Positive + fap_prop$FAP_Negative
fap_prop$Fraction_FAP_Pos <- fap_prop$FAP_Positive / fap_prop$Total

ggplot(fap_prop, aes(x = reorder(Cluster, -Fraction_FAP_Pos), y = Fraction_FAP_Pos)) +
  geom_bar(stat = "identity", fill = "firebrick") +
  coord_flip() +
  labs(title = "Proportion of FAP+ Cells per Cluster", y = "Fraction FAP+", x = "") +
  theme_minimal()

hist(FetchData(mpnst_combined, "FAP")$FAP, breaks = 100, main = "FAP Expression", xlab = "Normalized Expression")
abline(v = 0.5, col = "red", lty = 2)


# Marker gene sets for MPNST subtypes
g1_markers <- c("TWIST1", "SOX9", "SNAI2", "OTX2", "PAX3", "PAX6", "SMO")
g2_markers <- c("GAP43", "PLP1", "NGFR", "PTCH1", "WNT11")

# Adding module scores to the combined object
mpnst_combined <- AddModuleScore(mpnst_combined, features = list(g1_markers), name = "G1_score")
mpnst_combined <- AddModuleScore(mpnst_combined, features = list(g2_markers), name = "G2_score")

# tumor_state based on module scores
mpnst_combined$tumor_state <- ifelse(
  mpnst_combined$G1_score1 > mpnst_combined$G2_score1,
  "G1",
  "G2"
)

mpnst_combined$tumor_state <- factor(mpnst_combined$tumor_state)

#All cells
DimPlot(mpnst_combined, group.by = "tumor_state", label = TRUE) + ggtitle("Tumor State (G1 vs G2)")

# Define tumor clusters
tumor_clusters <- c("Tumor_Mesenchymal_like", "Tumor_Neurodevelopmental", "Cycling_Tumor", "Tumor_Schwann_like")

tumor_cells <- subset(mpnst_combined, idents = tumor_clusters)

tumor_cells$tumor_state <- factor(tumor_cells$tumor_state)

#UMAP colored by tumor_state (G1 vs G2)
DimPlot(tumor_cells, reduction = "umap", group.by = "tumor_state", 
             pt.size = 0.5) + 
  ggtitle("Tumor States (G1 vs G2) in Tumor Clusters") +
  theme(plot.title = element_text(hjust = 0.5))

#EMT module score
emt_genes <- strsplit(
  "ABI3BP,ACTA2,ADAM12,ANPEP,APLP1,AREG,BASP1,BDNF,BGN,BMP1,CADM1,CALD1,CALU,CAP2,CAPG,CD44,CD59,CDH11,CDH2,CDH6,COL11A1,COL12A1,COL16A1,COL1A1,COL1A2,COL3A1,COL4A1,COL4A2,COL5A1,COL5A2,COL5A3,COL6A2,COL6A3,COL7A1,COL8A2,COMP,COPA,CRLF1,CTGF,CTHRC1,CXCL1,CXCL12,CXCL6,CYR61,DAB2,DCN,DKK1,DPYSL3,DST,ECM1,ECM2,EDIL3,EFEMP2,ELN,EMP3,ENO2,FAP,FAS,FBLN1,FBLN2,FBLN5,FBN1,FBN2,FERMT2,FGF2,FLNA,FMOD,FN1,FOXC2,FSTL1,FSTL3,FUCA1,FZD8,GADD45A,GADD45B,GAS1,GEM,GJA1,GLIPR1,GLT25D1,GPC1,GPX7,GREM1,HTRA1,ID2,IGFBP2,IGFBP3,IGFBP4,IL15,IL32,IL6,IL8,INHBA,ITGA2,ITGA5,ITGAV,ITGB1,ITGB3,ITGB5,JUN,LAMA1,LAMA2,LAMA3,LAMC1,LAMC2,LEPRE1,LGALS1,LOX,LOXL1,LOXL2,LRP1,LRRC15,LUM,MAGEE1,MATN2,MATN3,MCM7,MEST,MFAP5,MGP,MMP1,MMP14,MMP2,MMP3,MSX1,MXRA5,MYL9,MYLK,NID2,NNMT,NOTCH2,NT5E,NTM,OXTR,PCOLCE,PCOLCE2,PDGFRB,PDLIM4,PFN2,PLAUR,PLOD1,PLOD2,PLOD3,PMEPA1,PMP22,POSTN,PPIB,PRRX1,PRSS2,PTHLH,PTX3,PVR,QSOX1,RGS4,RHOB,SAT1,SCG2,SDC1,SDC4,SERPINE1,SERPINE2,SERPINH1,SFRP1,SFRP4,SGCB,SGCD,SGCG,SLC6A8,SLIT2,SLIT3,SNAI2,SNTB1,SPARC,SPOCK1,SPP1,TAGLN,TFPI2,TGFB1,TGFBI,TGFBR3,TGM2,THBS1,THBS2,THY1,TIMP1,TIMP3,TNC,TNFAIP3,TNFRSF11B,TNFRSF12A,TPM1,TPM2,TPM4,VCAM1,VCAN,VEGFA,VEGFC,VIM,WIPF1,WNT5A",
  split = ","
)[[1]]

mpnst_combined<- AddModuleScore(
  object = mpnst_combined,
  features = list(emt_genes),
  name = "EMT_Score"
)

#Annotation
mpnst_combined$cell_type <- plyr::mapvalues(
  x = Idents(mpnst_combined),
  from = names(cluster_ids),
  to = cluster_ids
)

#Data
df <- FetchData(mpnst_combined, vars = c("cell_type", "EMT_Score1"))

ggplot(df, aes(x = cell_type, y = EMT_Score1)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 80, hjust = 1)) +
  labs(title = "EMT Score Across Cell Types", x = "Cell Type", y = "EMT Score")

tumor_clusters <- c(
  "Tumor_Mesenchymal_like",
  "Tumor_Schwann_like",
  "Tumor_Neurodevelopmental",
  "Cycling_Tumor"
)

df <- FetchData(mpnst_combined, vars = c("FAP", "EMT_Score1", "cell_type"))
df_tumor <- df %>% filter(cell_type %in% tumor_clusters)

# Plot EMT score across tumor cell types only
ggplot(df_tumor, aes(x = cell_type, y = EMT_Score1)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 80, hjust = 1)) +
  labs(title = "EMT Score Across Tumor Cell Types", x = "Tumor Cell Type", y = "EMT Score")

tumor_clusters <- c(
  "Tumor_Mesenchymal_like",
  "Tumor_Schwann_like",
  "Tumor_Neurodevelopmental",
  "Cycling_Tumor"
)

# Filter metadata for tumor clusters only
df <- FetchData(mpnst_combined, vars = c("FAP", "EMT_Score1", "cell_type"))
df_tumor <- df %>% filter(cell_type %in% tumor_clusters)

# Plot EMT score across tumor cell types only
ggplot(df_tumor, aes(x = cell_type, y = EMT_Score1)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 80, hjust = 1)) +
  labs(title = "EMT Score Across Tumor Cell Types", x = "Tumor Cell Type", y = "EMT Score")

#Spearman correlations of FAP vs EMT score within tumor clusters
cor_results_tumor <- df_tumor %>%
  group_by(cell_type) %>%
  summarise(
    spearman_rho = cor(FAP, EMT_Score1, method = "spearman", use = "complete.obs"),
    p_value = cor.test(FAP, EMT_Score1, method = "spearman")$p.value,
    .groups = "drop"
  )

print(cor_results_tumor)

Idents(mpnst_combined)
names(cluster_ids)[cluster_ids %in% tumor_clusters]
levels(Idents(mpnst_combined))
tumor_clusters <- c("Tumor_Mesenchymal_like", "Tumor_Neurodevelopmental", "Cycling_Tumor")

cells_tumor <- WhichCells(mpnst_combined, idents = tumor_clusters)

seurat_tumor <- subset(mpnst_combined, cells = cells_tumor)

FeaturePlot(
  object = seurat_tumor,
  features = "EMT_Score1",
  cols = c("white", "red"),
  reduction = "umap"
) + ggtitle("EMT Score across Tumor Clusters")

FeaturePlot(
  object = seurat_tumor,
  features = "FAP",
  cols = c("lightgrey", "red"),
  reduction = "umap"
) + ggtitle("FAP Expression")


# Add significance labels based on p-values
cor_results_tumor <- cor_results_tumor %>%
  mutate(
    signif = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

ggplot(cor_results_tumor, aes(x = cell_type, y = spearman_rho, fill = spearman_rho)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = signif), vjust = -0.2, size = 6) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(
    title = "FAP and EMT Score across Tumor Clusters",
    x = "Tumor Cluster",
    y = "Spearman Correlation (rho)"
  ) +
  coord_cartesian(ylim = c(0, 0.4)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Differentiation trajector from neural-crest to Schwann-like
lineage_clusters <- c("Tumor_Neurodevelopmental", "Tumor_Schwann_like")
cells_lineage   <- WhichCells(mpnst_combined, idents = lineage_clusters)
lineage_subset  <- subset(mpnst_combined, cells = cells_lineage)

#Umap and Cluster labels
umap_emb <- Embeddings(lineage_subset, "umap")
clust     <- Idents(lineage_subset)

#Slingshoth with start at neural-crest like and end at Schwann-like
ss <- slingshot(umap_emb,
                clusterLabels = clust,
                start.clus     = "Tumor_Neurodevelopmental",
                end.clus       = "Tumor_Schwann_like")

#Pseudotime
pseudo_nc2scp <- slingPseudotime(ss)[,1]

lineage_subset$pseudotime_nc2scp <- pseudo_nc2scp
FeaturePlot(lineage_subset, "pseudotime_nc2scp", reduction = "umap") +
  scale_color_viridis_c() +
  ggtitle("Neural crest to Schwann-like/precursor")

#FAP expression along that trajectory
df <- data.frame(
  pseudotime = lineage_subset$pseudotime_nc2scp,
  FAP        = FetchData(lineage_subset, "FAP")[,1]
)
df <- na.omit(df)

ggplot(df, aes(pseudotime, FAP)) +
  geom_point(alpha=0.5) +
  geom_smooth(method="loess") +
  theme_minimal() +
  labs(title="FAP vs Neural crest to Schwann",
       x="Pseudotime", y="FAP expression")

cor.test(df$pseudotime, df$FAP, method="spearman")


df <- data.frame(
  pseudotime = lineage_subset$pseudotime_nc2scp,
  FAP        = FetchData(lineage_subset, "FAP")[,1]
)
df <- na.omit(df)

#Spearman rho
rho_val <- cor.test(df$pseudotime, df$FAP, method = "spearman")$estimate

ggplot(df, aes(x = pseudotime, y = FAP)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_smooth(method = "loess", se = TRUE) +
  theme_minimal() +
  labs(
    title    = sprintf("FAP and differentiation trajectory (ρ=%.2f)", rho_val),
    x        = "Pseudotime",
    y        = "FAP expression"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

#Plot FAP expression
FeaturePlot(lineage_subset, "FAP", reduction = "umap") +
  scale_color_viridis_c() +
  ggtitle("FAP: Neural crest to Schwann-like/precursor")

#Plot SMO expression
FeaturePlot(lineage_subset, "SMO", reduction = "umap") +
  scale_color_viridis_c() +
  ggtitle("SMO: Neural crest to Schwann-like/precursor")


sessionInfo()
