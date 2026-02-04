# ============================================================
# Title: Tissue-Specific Tau Heatmap (Ribo-Seq & Proteogenomics)
# Author: Brendan Miller
# Description:
#   Filters high-confidence microproteins with evidence support
#   (Ribo-seq / SEER / SAM), applies tau > 0.75 tissue-specific
#   threshold, and visualizes log2CPM profiles as Muskified
#   publication-ready heatmaps.
# ============================================================

# === Libraries ===
library(dplyr)
library(readr)
library(pheatmap)
library(stringr)
library(ggplot2)
library(Cairo)

# ============================================================
# === Load Data ===
# ============================================================
df <- read_csv(
  "/Users/brendanmiller/Library/CloudStorage/Box-Box/panda/dictionaries/sequence_hits_dictionary.csv"
)

# ============================================================
# === Evidence & Filtering Logic ===
# ============================================================
df <- df %>%
  mutate(
    # 1. Base Evidence Flags
    has_riboseq = !is.na(riboseq_source_count) & riboseq_source_count > 0,
    has_seer    = !is.na(proteogenomics) & proteogenomics == "SEER",
    has_sam     = str_detect(coalesce(shortstop_predictions, ""), "SAM"),
    is_swiss    = coalesce(primary_smorf_type == "Swiss-Prot-MP", FALSE),
    
    # 2. Shared Type Restriction
    # Must NOT be oORF, isoORF, Unknown, or Empty
    pass_type_check = !(primary_smorf_type %in% c("oORF", "isoORF", "Unknown", "") | is.na(primary_smorf_type)),
    
    # 3. Tightened CPM Restriction (Requires at least one tissue to be > 5, NOT just NA)
    # We use coalesce(..., 0) so that NAs are treated as 0 (failing the > 5 test)
    pass_cpm_check = coalesce(muscle_riboseq_cpm_per_kb, 0) > 5 | 
      coalesce(heart_riboseq_cpm_per_kb, 0)  > 5 | 
      coalesce(fat_riboseq_cpm_per_kb, 0)    > 5 | 
      coalesce(liver_riboseq_cpm_per_kb, 0)   > 5
  ) %>%
  filter(
    is_swiss | (
      coalesce(protein_class == "Microprotein", FALSE) & (
        # Door 1: It's SEER (Independent of types/CPM)
        has_seer | 
          
          # Door 2: It's Ribo-seq AND has a valid type
          (has_riboseq & pass_type_check) | 
          
          # Door 3: It's SAM AND has a valid type AND passes CPM
          (has_sam & pass_cpm_check)
      )
    )
  ) %>%
  arrange(desc(is_swiss)) %>%
  distinct(sequence, .keep_all = TRUE)


# ============================================================
# === Tau Filtering & Matrix Prep ===
# ============================================================
tau_expression <- df %>%
  filter(tau_score > 0.9)
dim(tau_expression)
expression_matrix <- tau_expression %>%
  dplyr::select(starts_with("log2CPM_")) %>%
  as.data.frame()

# Clean column names (remove "log2CPM_" prefix)
colnames(expression_matrix) <- gsub("^log2CPM_", "", colnames(expression_matrix))

# ============================================================
# === Muskified Heatmap — Global Expression ===
# ============================================================
heatmap_colors <- colorRampPalette(c("white","#e8f4ff", "#b22222"))(100)

# heatmap1 <- pheatmap(
#   expression_matrix,
#   scale = "row",
#   cluster_rows = TRUE,
#   cluster_cols = TRUE,
#   show_rownames = FALSE,
#   show_colnames = TRUE,
#   color = heatmap_colors,
#   clustering_distance_rows = "euclidean",
#   clustering_distance_cols = "euclidean",
#   clustering_method = "complete",
#   border_color = NA,
#   fontsize = 12,
#   fontsize_col = 14,
#   legend = TRUE,
#   main = "Tau > 0.75 | Rhythymic Transcription of smORFs",
#   angle_col = 45
# )
# heatmap1
# # === Save Heatmap ===
# ggsave(
#   filename = "/Users/brendanmiller/Library/CloudStorage/Box-Box/brendan_alan_shared_folder/Panda_paper/figures/tau_heatmap.pdf",
#   plot = heatmap1,
#   width = 6, height = 4, device = cairo_pdf)
# 
# grid::grid.newpage()
# grid::grid.draw(heatmap1$gtable)
# dev.off()

# ============================================================
# === Muskified Heatmap — Clustered (k-means = 11) ===
# ============================================================

# Set seed for reproducible k-means clustering
set.seed(42)

heatmap2 <- pheatmap(
  expression_matrix,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  kmeans_k = 11,
  show_rownames = TRUE,
  show_colnames = TRUE,
  color = heatmap_colors,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  #border_color = NA,
  fontsize = 12,
  fontsize_col = 14,
  legend = TRUE,
  main = "Tau > 0.9 | k-means Clusters (k=11)",
  angle_col = 45
)

# === Save Clustered Heatmap ===
ggsave(
  filename = "/Users/brendanmiller/Library/CloudStorage/Box-Box/brendan_alan_shared_folder/Panda_paper/figures/figure2/panels/tau_heatmap.pdf",
  plot = heatmap2,
  width = 9, height = 4, device = cairo_pdf)

# === Save pheatmap object for future sessions ===
saveRDS(
  heatmap2,
  "/Users/brendanmiller/Library/CloudStorage/Box-Box/brendan_alan_shared_folder/Panda_paper/supplemental/tau_heatmap_object.rds"
)
cat("Saved heatmap object to tau_heatmap_object.rds\n")
cat("To reload in future: heatmap2 <- readRDS('supplemental/tau_heatmap_object.rds')\n")

# ============================================================
# === Export Supplementary CSV with Tau Clusters ===
# ============================================================

# Extract k-means cluster assignments from pheatmap object
kmeans_clusters <- heatmap2$kmeans$cluster

# Create supplementary table with requested columns
tau_supplementary <- tau_expression %>%
  mutate(tau_cluster = kmeans_clusters) %>%
  dplyr::select(
    CLICK_UCSC,
    sequence,
    primary_smorf_type,
    tau_cluster,
    gene_ids
  ) %>%
  arrange(tau_cluster)

# Save to CSV
write_csv(
  tau_supplementary,
  "/Users/brendanmiller/Library/CloudStorage/Box-Box/brendan_alan_shared_folder/Panda_paper/supplemental/S1_figure2_tau_output.csv"
)

# Print summary
cat("Supplementary CSV saved with", nrow(tau_supplementary), "rows\n")
cat("Cluster distribution:\n")
print(table(tau_supplementary$tau_cluster))