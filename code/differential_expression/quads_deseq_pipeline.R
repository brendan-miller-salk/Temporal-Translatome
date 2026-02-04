#!/usr/bin/env Rscript
#===============================================================================
# DESeq2: Differential Expression Analysis for PANDA Mouse Atlas
#===============================================================================
# This script documents the DESeq2 pipeline used to identify differentially
# expressed genes and microproteins in mouse quadriceps muscle.
#
# Experimental Design:
#   - Tissue: Quadriceps muscle
#   - Factors: Sex (M/F) × Exercise (Sed/RW) × Time (ZT4/ZT16)
#   - Full factorial design with all interactions
#
# Tool: DESeq2 v1.38+
# Reference: Love MI, Huber W, Anders S. (2014) Genome Biology 15:550.
#            PMID: 25516281
#===============================================================================

# ==============================================================================
# SETUP
# ==============================================================================

#-------------------------------------------------------------------------------
# Load Required Libraries
#-------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(biomaRt)
  library(BiocParallel)
  library(data.table)
})

#-------------------------------------------------------------------------------
# Configure Parallel Processing
#-------------------------------------------------------------------------------

register(MulticoreParam(workers = 8))

# ==============================================================================
# DATA LOADING
# ==============================================================================

#-------------------------------------------------------------------------------
# Load Count Matrix
#-------------------------------------------------------------------------------
# Input: Gene-by-sample count matrix from featureCounts
# Rows: Gene/transcript IDs (ENSMUSG* for canonical, MSTRG* for novel)
# Columns: Sample names

counts <- fread("combined_counts.csv", data.table = FALSE)
rownames(counts) <- counts[[1]]
counts <- counts[, -1]

# Clean sample names (remove alignment suffixes)
colnames(counts) <- gsub("Aligned\\.sortedByCoord\\.out\\.bam", "", colnames(counts))
counts <- counts[, !duplicated(colnames(counts))]

#-------------------------------------------------------------------------------
# Load and Process Sample Metadata
#-------------------------------------------------------------------------------
# Required columns: FileName, Sex, Group (exercise), Time, Tissue

metadata <- read.csv("sample_metadata.csv") %>%
  filter(Method == "RNASeq") %>%
  mutate(FileName = gsub("_R[12]_001\\.fastq\\.gz", "", FileName)) %>%
  distinct(FileName, .keep_all = TRUE)

# Align metadata with count matrix columns
metadata <- metadata[match(colnames(counts), metadata$FileName), ] %>%
  filter(complete.cases(FileName))

# Set reference levels for each factor
metadata <- metadata %>%
  mutate(
    Sex   = relevel(factor(Sex),   ref = "M"),      # Male as reference
    Group = relevel(factor(Group), ref = "Sed"),    # Sedentary as reference
    Time  = relevel(factor(Time),  ref = "ZT16")    # ZT16 as reference
  )

# Final alignment
counts <- counts[, intersect(metadata$FileName, colnames(counts))]
rownames(metadata) <- colnames(counts)

# ==============================================================================
# DESEQ2 ANALYSIS
# ==============================================================================

#-------------------------------------------------------------------------------
# Create DESeqDataSet
#-------------------------------------------------------------------------------
# Full factorial design: Sex × Exercise × Time
# This models main effects and all 2-way and 3-way interactions

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = metadata,
  design    = ~ Sex * Group * Time
)

#-------------------------------------------------------------------------------
# Filter Low-Count Genes
#-------------------------------------------------------------------------------
# Remove genes with fewer than 5 total counts across all samples

dds <- dds[rowSums(counts(dds)) >= 5, ]
message(sprintf("Genes after filtering: %d", nrow(dds)))

#-------------------------------------------------------------------------------
# Run DESeq2
#-------------------------------------------------------------------------------
# fitType = "local": Use local regression for dispersion estimation
# sfType = "poscounts": Handle genes with zero counts in some samples

dds <- DESeq(
  dds,
  fitType  = "local",
  sfType   = "poscounts",
  parallel = TRUE
)

# Save fitted object for reproducibility
saveRDS(dds, "dds_fitted.rds")

# ==============================================================================
# RESULTS EXTRACTION
# ==============================================================================

#-------------------------------------------------------------------------------
# Extract All Contrasts of Interest
#-------------------------------------------------------------------------------

res_list <- list(
  # Main effects
  Exercise_main = results(dds, name = "Group_RW_vs_Sed"),
  Time_main     = results(dds, name = "Time_ZT4_vs_ZT16"),
  Sex_main      = results(dds, name = "Sex_F_vs_M"),
  
  # Two-way interactions

Exercise_x_Time = results(dds, name = "GroupRW.TimeZT4"),
  Sex_x_Exercise  = results(dds, name = "SexF.GroupRW"),
  Sex_x_Time      = results(dds, name = "SexF.TimeZT4"),
  
  # Three-way interaction
  Sex_x_Exercise_x_Time = results(dds, name = "SexF.GroupRW.TimeZT4"),
  
  # Simple effects (exercise effect at specific time points)
  Exercise_at_ZT4  = results(dds, contrast = list(c("Group_RW_vs_Sed", "GroupRW.TimeZT4"))),
  Exercise_at_ZT16 = results(dds, name = "Group_RW_vs_Sed")
)

#-------------------------------------------------------------------------------
# Combine Results into Single Table
#-------------------------------------------------------------------------------

results_combined <- lapply(names(res_list), function(name) {
  res_list[[name]] %>%
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    dplyr::select(gene_id, baseMean, log2FoldChange, pvalue, padj) %>%
    dplyr::rename_with(~ paste0(name, "_", .), -1)
})

combined_df <- Reduce(function(x, y) full_join(x, y, by = "gene_id"), results_combined)

#-------------------------------------------------------------------------------
# Filter for Significant Genes
#-------------------------------------------------------------------------------
# Keep genes significant (padj < 0.05) in any contrast

sig_results <- combined_df %>%
  filter(if_any(ends_with("padj"), ~ .x < 0.05))

message(sprintf("Significant genes: %d", nrow(sig_results)))

# ==============================================================================
# GENE ANNOTATION
# ==============================================================================

#-------------------------------------------------------------------------------
# Add Gene Symbols via biomaRt
#-------------------------------------------------------------------------------

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Extract Ensembl IDs (strip version numbers)
sig_results <- sig_results %>%
  mutate(ensembl_id = ifelse(
    grepl("^ENSMUSG", gene_id),
    sub("\\..*", "", gene_id),
    NA_character_
  ))

# Query biomaRt for gene symbols
annot <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters    = "ensembl_gene_id",
  values     = na.omit(unique(sig_results$ensembl_id)),
  mart       = ensembl
) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

# Merge annotations
sig_results <- sig_results %>%
  left_join(annot, by = c("ensembl_id" = "ensembl_gene_id")) %>%
  relocate(mgi_symbol, .after = gene_id) %>%
  dplyr::select(-ensembl_id)

# ==============================================================================
# EXPORT RESULTS
# ==============================================================================

write.csv(sig_results, "differential_expression_results.csv", row.names = FALSE)

message("Analysis complete. Results saved to differential_expression_results.csv")

# ==============================================================================
# VISUALIZATION FUNCTIONS
# ==============================================================================

#-------------------------------------------------------------------------------
# Bar Plot with Error Bars
#-------------------------------------------------------------------------------
# Creates publication-ready bar plots showing normalized counts by condition

plot_counts_bar <- function(gene_id, dds_obj, intgroup = c("Group", "Time", "Sex")) {
  
  stopifnot(gene_id %in% rownames(dds_obj))
  
  # Extract normalized counts
  df <- plotCounts(dds_obj, gene = gene_id, intgroup = intgroup, returnData = TRUE)
  
  # Set factor levels
  df <- df %>%
    mutate(
      Group = factor(Group, levels = c("Sed", "RW")),
      Time  = factor(Time,  levels = c("ZT16", "ZT4")),
      Sex   = factor(Sex,   levels = c("M", "F"))
    )
  
  # Calculate summary statistics
  summary_df <- df %>%
    group_by(Group, Time, Sex) %>%
    summarise(
      mean = mean(count),
      sem  = sd(count) / sqrt(n()),
      .groups = "drop"
    )
  
  # Create plot
  ggplot(summary_df, aes(x = Sex, y = mean, fill = Group)) +
    geom_col(
      position = position_dodge(width = 0.6),
      width = 0.5, alpha = 0.7,
      color = "black", linewidth = 0.5
    ) +
    geom_errorbar(
      aes(ymin = mean - sem, ymax = mean + sem),
      position = position_dodge(width = 0.6),
      width = 0.15, linewidth = 0.5
    ) +
    geom_point(
      data = df,
      aes(x = Sex, y = count, fill = Group),
      position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.6),
      size = 3, shape = 21, stroke = 0.3
    ) +
    facet_wrap(~ Time) +
    scale_fill_manual(
      values = c("Sed" = "#808080", "RW" = "#C41E3A"),
      labels = c("Sed" = "Sedentary", "RW" = "Running Wheel")
    ) +
    labs(
      title = gene_id,
      x = NULL,
      y = "Normalized Counts",
      fill = "Exercise"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA),
      panel.grid = element_blank(),
      legend.position = "bottom",
      strip.text = element_text(face = "bold")
    )
}

# ==============================================================================
# EXAMPLE USAGE
# ==============================================================================

# # Load pre-computed DESeq object
# dds <- readRDS("dds_fitted.rds")
#
# # Plot a gene of interest
# p <- plot_counts_bar("ENSMUSG00000026051", dds)
# ggsave("example_gene_plot.pdf", p, width = 8, height = 5)
#
# # Get results for a specific contrast
# exercise_results <- results(dds, name = "Group_RW_vs_Sed")
# exercise_sig <- subset(exercise_results, padj < 0.05)
