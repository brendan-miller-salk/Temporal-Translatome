# ============================================================
# 🧬 Microprotein smORF Type Visualization & Tissue Overlap
# Project: Panda Atlas Analysis
# Author: Brendan Miller
# Date: 2025-10
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(reticulate)
})

# Optional venn (install.packages("ggvenn") if missing)
has_ggvenn <- requireNamespace("ggvenn", quietly = TRUE)

# ---------------------------
# Config (edit for your env)
# ---------------------------
cfg <- list(
  csv_path = "/Users/brendanmiller/Library/CloudStorage/Box-Box/panda/dictionaries/sequence_hits_dictionary.csv",
  save_dir = "//Users/brendanmiller/Library/CloudStorage/Box-Box/brendan_alan_shared_folder/Panda_paper/figures/figure1/panels",
  id_col   = "sequence"  # unique microprotein identifier
)

# ---------------------------
# Utils
# ---------------------------
log_info <- function(...) cat(sprintf(paste0(..., "\n")))
ensure_cols <- function(df, cols) {
  missing <- setdiff(cols, colnames(df))
  if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

# ---------------------------
# Theme & Palette
# ---------------------------
colors_fill <- c(
  "Swiss-Prot-MP" = "#9EC9E2",  # soft blue
  "Noncanonical"  = "#FDE28D"   # pale gold
)

custom_theme <- theme_minimal(base_family = "roboto") +
  theme(
    plot.title       = element_text(hjust = 0.1, size = 30),
    axis.text.x      = element_text(hjust = 1, size = 25, color = "black"),
    axis.text.y      = element_text(size = 25, color = "black"),
    axis.title       = element_text(size = 25, margin = margin(r = 10)),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "none"
  )

# ---------------------------
# Load & Filter
# ---------------------------
load_data <- function(path) {
  if (!file.exists(path)) stop("CSV not found: ", path)
  readr::read_csv(path, show_col_types = FALSE) %>% as_tibble()
}

filter_df <- function(df) {
  # Required columns for filtering
  ensure_cols(df, c(
    "protein_class","primary_smorf_type","riboseq_source_count","proteogenomics",
    "shortstop_predictions","muscle_riboseq_cpm_per_kb","heart_riboseq_cpm_per_kb",
    "fat_riboseq_cpm_per_kb","liver_riboseq_cpm_per_kb"
  ))
  
  df %>%
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
  
}

# ---------------------------
# Tissue Logic (CPM gating for non-riboseq)
# ---------------------------
add_tissue_flags <- function(df) {
  ensure_cols(df, c("has_riboseq","has_sam",
                    "has_liver_riboseq","has_fat_riboseq","has_heart_riboseq","has_muscle_riboseq",
                    "liver_riboseq_cpm_per_kb","fat_riboseq_cpm_per_kb",
                    "heart_riboseq_cpm_per_kb","muscle_riboseq_cpm_per_kb"))
  
  df %>%
    mutate(across(starts_with("has_"), ~ as.logical(.x))) %>%
    mutate(
      in_liver  = has_liver_riboseq  | (!has_riboseq & has_sam & liver_riboseq_cpm_per_kb  > 5),
      in_fat    = has_fat_riboseq    | (!has_riboseq & has_sam & fat_riboseq_cpm_per_kb    > 5),
      in_heart  = has_heart_riboseq  | (!has_riboseq & has_sam & heart_riboseq_cpm_per_kb  > 5),
      in_muscle = has_muscle_riboseq | (!has_riboseq & has_sam & muscle_riboseq_cpm_per_kb > 5)
    ) %>%
    rowwise() %>%
    mutate(tissues_detected = sum(c_across(c(in_liver,in_fat,in_heart,in_muscle)), na.rm = TRUE)) %>%
    ungroup()
}

# ---------------------------
# Summaries
# ---------------------------
summaries <- function(df) {
  # Evidence breakdown table
  evidence_summary <- df %>%
    mutate(Evidence = case_when(
      primary_smorf_type == "Swiss-Prot-MP" ~ "Swiss-Prot-MP",
      has_riboseq & has_sam & has_seer ~ "Ribo+SAM+SEER",
      has_riboseq & has_sam ~ "Ribo+SAM",
      has_riboseq & has_seer ~ "Ribo+SEER",
      has_sam & has_seer ~ "SAM+SEER",
      has_riboseq ~ "Ribo-only",
      has_sam ~ "SAM-only",
      has_seer ~ "SEER-only",
      TRUE ~ "Other"
    )) %>%
    count(Evidence, name = "count") %>%
    arrange(desc(count))
  
  canonical_summary <- df %>%
    mutate(category = if_else(primary_smorf_type == "Swiss-Prot-MP", "Swiss-Prot-MP", "Noncanonical")) %>%
    count(category, name = "count")
  
  tissue_summary <- df %>%
    count(tissues_detected, name = "microprotein_count") %>%
    mutate(category = case_when(
      tissues_detected == 0 ~ "None",
      tissues_detected == 1 ~ "Tissue-specific",
      tissues_detected == 2 ~ "Shared in 2 tissues",
      tissues_detected == 3 ~ "Shared in 3 tissues",
      tissues_detected >= 4 ~ "Shared in all 4 tissues"
    ))
  
  tissue_counts <- df %>%
    summarise(
      liver  = sum(in_liver,  na.rm = TRUE),
      fat    = sum(in_fat,    na.rm = TRUE),
      heart  = sum(in_heart,  na.rm = TRUE),
      muscle = sum(in_muscle, na.rm = TRUE)
    ) %>%
    pivot_longer(everything(), names_to = "tissue", values_to = "count") %>%
    arrange(desc(count))
  
  list(
    evidence = evidence_summary,
    canonical = canonical_summary,
    tissue_overlap = tissue_summary,
    tissue_counts = tissue_counts
  )
}

# ---------------------------
# Plots
# ---------------------------
# ---------------------------
# Enhanced smORF Type Plot
# ---------------------------
plot_smorf_bar <- function(df) {
  library(ggtext)
  
  smorf_counts <- df %>%
    count(primary_smorf_type, name = "count") %>%
    mutate(
      bar_fill = if_else(primary_smorf_type == "Swiss-Prot-MP", "Swiss-Prot-MP", "Noncanonical")
    ) %>%
    arrange(desc(count))
  
  # --- Metallic Gradient Palette ---
  fill_palette <- c("Swiss-Prot-MP" = "#A9C9E2", "Noncanonical" = "#FAD97F")
  
  # --- Plot ---
  ggplot(smorf_counts, aes(x = count, y = reorder(primary_smorf_type, count), fill = bar_fill)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.8) +
    scale_fill_manual(values = fill_palette, na.value = "gray80") +
    geom_text(
      aes(label = scales::comma(count)),
      hjust = -0.1,
      size = 7,
      family = "roboto",
      color = "black"
    ) +
    labs(
      title = "Distribution of <span style='color:#34607b;'>Swiss-Prot</span> vs <span style='color:#a28331;'>Noncanonical</span> smORFs",
      x = "Count",
      y = NULL
    ) +
    custom_theme +
    theme(
      plot.title       = element_markdown(size = 28, hjust = 0.05),
      axis.text.y      = element_text(size = 24, color = "black"),
      axis.text.x      = element_text(size = 22, color = "black"),
      axis.title.x     = element_text(size = 22),
      panel.background = element_rect(fill = "#FBFBFB", color = NA),
      panel.grid.major.x = element_line(linetype = "dashed", color = "gray80"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 30, 10, 20)
    ) +
    coord_cartesian(expand = TRUE, clip = "off") +
    guides(fill = "none")
}

plot_venn <- function(df, id_col = "sequence") {
  if (!has_ggvenn) {
    log_info("ggvenn not installed; skipping Venn. install.packages('ggvenn') to enable.")
    return(NULL)
  }
  if (!id_col %in% colnames(df)) stop("Identifier column not found: ", id_col)
  
  venn_list <- list(
    Liver  = df %>% filter(in_liver)  %>% pull(all_of(id_col)),
    Fat    = df %>% filter(in_fat)    %>% pull(all_of(id_col)),
    Heart  = df %>% filter(in_heart)  %>% pull(all_of(id_col)),
    Muscle = df %>% filter(in_muscle) %>% pull(all_of(id_col))
  )
  
  ggvenn::ggvenn(
    venn_list,
    fill_color   = c("#FBE7A1", "#F6C453", "#EBA937", "#C87E1C"),
    stroke_color = "black",
    stroke_size  = 0.5,
    set_name_size = 8,
    text_size    = 5
  ) +
    labs(title = "Shared Microproteins Across Tissues") +
    theme_minimal(base_family = "roboto") +
    theme(
      plot.title       = element_text(hjust = 0.5, size = 30),
      axis.text.x      = element_blank(),
      axis.text.y      = element_blank(),
      axis.title       = element_blank(),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position  = "none"
    )
}

# ---------------------------
# Save helpers
# ---------------------------
save_plot <- function(p, path_pdf, path_png = NULL, w = 10, h = 8) {
  dir.create(dirname(path_pdf), recursive = TRUE, showWarnings = FALSE)
  ggsave(path_pdf, plot = p, width = w, height = h, units = "in", dpi = 300, device = cairo_pdf)
  if (!is.null(path_png)) ggsave(path_png, plot = p, width = w, height = h, units = "in", dpi = 300)
}

save_supplementary <- function(df_filtered, dirpath) {
  dir.create(dirpath, recursive = TRUE, showWarnings = FALSE)
  
  # Create single supplementary file with specified columns
  supplementary <- df_filtered %>%
    dplyr::select(
      UCSC_link = CLICK_UCSC,
      sequence,
      gene_id = gene_ids,
      gene_name = gene_body,
      smORF_type = primary_smorf_type,
      has_liver_riboseq,
      has_fat_riboseq,
      has_heart_riboseq,
      has_muscle_riboseq,
      shortstop_predictions,
      shortstop_confidences,
      has_seer,
      muscle_riboseq_cpm_per_kb,
      heart_riboseq_cpm_per_kb,
      fat_riboseq_cpm_per_kb,
      liver_riboseq_cpm_per_kb
    )
  
  out_path <- "/Users/brendanmiller/Library/CloudStorage/Box-Box/brendan_alan_shared_folder/Panda_paper/supplemental/S1_figure1_output.csv"
  readr::write_csv(supplementary, out_path)
  log_info("Saved supplementary file to ", out_path, " with ", nrow(supplementary), " microproteins.")
}

# ---------------------------
# Supplementary Figure Functions
# ---------------------------

# SFig1: Amino Acid Frequency Bar Graph
plot_aa_frequency <- function(df) {
  # 20 standard amino acids only
  standard_aa <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
                   "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  
  # Add start_codon_type column and group label
  df_annotated <- df %>%
    mutate(
      start_codon_type = case_when(
        str_detect(coalesce(shortstop_predictions, ""), "SAM") ~ "SAM",
        str_detect(coalesce(shortstop_predictions, ""), "PRISM") ~ "PRISM",
        TRUE ~ "NonATG"
      ),
      # Create 4 categories: Swiss-Prot as one, Non-Swiss-Prot split by start codon
      category = if_else(
        primary_smorf_type == "Swiss-Prot-MP", 
        "Swiss-Prot", 
        paste0("Non-Swiss-Prot (", start_codon_type, ")")
      )
    )
  
  # Calculate AA frequencies - using strsplit and unnest_longer
  aa_freq <- df_annotated %>%
    filter(!is.na(sequence), sequence != "") %>%
    select(category, sequence) %>%
    # Split each sequence into individual characters
    mutate(amino_acids = strsplit(sequence, "")) %>%
    select(-sequence) %>%
    unnest_longer(amino_acids) %>%
    filter(amino_acids %in% standard_aa) %>%
    count(category, amino_acids) %>%
    group_by(category) %>%
    mutate(freq = n / sum(n)) %>%
    ungroup()
  
  ggplot(aa_freq, aes(x = amino_acids, y = freq, fill = category)) +
    geom_col(position = position_dodge(width = 0.9), color = "black", linewidth = 0.3) +
    scale_fill_manual(
      values = c(
        "Swiss-Prot" = "#A9C9E2",
        "Non-Swiss-Prot (SAM)" = "#E6A500", 
        "Non-Swiss-Prot (PRISM)" = "#F0C040", 
        "Non-Swiss-Prot (NonATG)" = "#FAD97F"
      ),
      name = "Category"
    ) +
    labs(
      title = "Amino Acid Frequency",
      x = "Amino Acid",
      y = "Frequency"
    ) +
    custom_theme +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18),
      legend.position = "bottom",
      legend.text = element_text(size = 12)
    )
}

# SFig2: Amino Acid Length Histogram
plot_aa_length <- function(df) {
  df_annotated <- df %>%
    filter(!is.na(sequence), sequence != "") %>%
    mutate(
      aa_length = nchar(sequence),
      start_codon_type = case_when(
        str_detect(coalesce(shortstop_predictions, ""), "SAM") ~ "SAM",
        str_detect(coalesce(shortstop_predictions, ""), "PRISM") ~ "PRISM",
        TRUE ~ "NonATG"
      ),
      swiss_status = if_else(primary_smorf_type == "Swiss-Prot-MP", "Swiss-Prot", "Non-Swiss-Prot"),
      # For Swiss-Prot, ignore start codon type; for Non-Swiss-Prot, use it
      category = if_else(swiss_status == "Swiss-Prot", "Swiss-Prot", paste0("Non-Swiss-Prot (", start_codon_type, ")"))
    )
  
  ggplot(df_annotated, aes(x = aa_length, fill = category)) +
    geom_histogram(position = "identity", alpha = 0.6, bins = 50, color = "black", linewidth = 0.3) +
    facet_wrap(~ swiss_status, ncol = 1, scales = "free_y") +
    scale_fill_manual(
      values = c(
        "Swiss-Prot" = "#A9C9E2",
        "Non-Swiss-Prot (SAM)" = "#E6A500", 
        "Non-Swiss-Prot (PRISM)" = "#F0C040", 
        "Non-Swiss-Prot (NonATG)" = "#FAD97F"
      ),
      name = "Category"
    ) +
    labs(
      title = "Amino Acid Length Distribution",
      x = "Amino Acid Length",
      y = "Count"
    ) +
    custom_theme +
    theme(legend.position = "bottom", legend.text = element_text(size = 12))
}

# SFig3: CPM Distribution by smORF Type
plot_cpm_by_smorf <- function(df) {
  cpm_data <- df %>%
    filter(primary_smorf_type != "Swiss-Prot-MP") %>%
    select(primary_smorf_type, 
           muscle_riboseq_cpm_per_kb, 
           heart_riboseq_cpm_per_kb, 
           fat_riboseq_cpm_per_kb, 
           liver_riboseq_cpm_per_kb) %>%
    pivot_longer(cols = ends_with("_cpm_per_kb"), 
                 names_to = "tissue", 
                 values_to = "cpm") %>%
    filter(!is.na(cpm), cpm > 0) %>%
    mutate(
      tissue = str_remove(tissue, "_riboseq_cpm_per_kb"),
      tissue = str_to_title(tissue),
      log_cpm = log10(cpm + 1)
    )
  
  ggplot(cpm_data, aes(x = reorder(primary_smorf_type, log_cpm, median), y = log_cpm, fill = primary_smorf_type)) +
    geom_boxplot(outlier.size = 0.5, color = "black", linewidth = 0.5) +
    facet_wrap(~ tissue, ncol = 4) +
    scale_fill_brewer(palette = "Set3") +
    labs(
      title = "Ribo-seq CPM Distribution by smORF Type",
      x = NULL,
      y = "log10(RPKM)"
    ) +
    custom_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      legend.position = "none"
    )
}

# ---------------------------
# Main
# ---------------------------
main <- function(cfg) {
  log_info("\n=== LOADING CSV ===")
  df_raw <- load_data(cfg$csv_path)
  ensure_cols(df_raw, c(cfg$id_col, "protein_class","primary_smorf_type"))
  
  log_info("\n=== FILTERING ===")
  df_filt <- filter_df(df_raw)
  
  log_info("\nEntries (filtered): ", nrow(df_filt))
  log_info("Unique sequences: ", dplyr::n_distinct(df_filt[[cfg$id_col]]))
  
  log_info("\n=== TISSUE FLAGS (with ShortStop CPM gating) ===")
  df_tissue <- add_tissue_flags(df_filt)
  
  # --- Summaries ---
  log_info("\n=== SUMMARIES ===")
  sums <- summaries(df_tissue)
  print(sums$tissue_overlap)
  print(sums$tissue_counts)
  
  # --- Bar Plot ---
  log_info("\n=== PLOTTING smORF DISTRIBUTION ===")
  p_bar <- plot_smorf_bar(df_tissue)
  print(p_bar)
  save_plot(
    p_bar,
    path_pdf = file.path(cfg$save_dir, "microprotein_smorf_distribution.pdf"),
    path_png = file.path(cfg$save_dir, "microprotein_smorf_distribution.png")
  )
  log_info("Saved smORF distribution plots.")
  
  # --- Venn Plot (optional) ---
  log_info("\n=== VENN DIAGRAM ===")
  p_venn <- plot_venn(df_tissue, id_col = cfg$id_col)
  if (!is.null(p_venn)) {
    print(p_venn)
    save_plot(
      p_venn,
      path_pdf = file.path(cfg$save_dir, "microprotein_tissue_venn.pdf"),
      path_png = file.path(cfg$save_dir, "microprotein_tissue_venn.png"),
      w = 8, h = 7
    )
    log_info("Saved Venn diagram.")
  }
  
  # --- Save supplementary ---
  log_info("\n=== SAVING SUPPLEMENTARY TABLE ===")
  save_supplementary(df_tissue, file.path(cfg$save_dir))
  
  # --- Supplementary Figures ---
  log_info("\n=== CREATING SUPPLEMENTARY FIGURES ===")
  
  # SFig1: Amino Acid Frequency
  p_aa_freq <- plot_aa_frequency(df_tissue)
  save_plot(
    p_aa_freq,
    path_pdf = file.path(cfg$save_dir, "SFig1_aa_frequency.pdf"),
    w = 13, h = 6
  )
  log_info("Saved SFig1: Amino acid frequency.")
  
  # SFig2: Amino Acid Length Distribution
  p_aa_length <- plot_aa_length(df_tissue)
  save_plot(
    p_aa_length,
    path_pdf = file.path(cfg$save_dir, "SFig2_aa_length.pdf"),
    w = 13, h = 6
  )
  log_info("Saved SFig2: Amino acid length distribution.")
  
  # SFig3: CPM Distribution by smORF Type
  p_cpm_smorf <- plot_cpm_by_smorf(df_tissue)
  save_plot(
    p_cpm_smorf,
    path_pdf = file.path(cfg$save_dir, "SFig3_cpm_smorf.pdf"),
    w = 13, h = 6
  )
  log_info("Saved SFig3: CPM by smORF type.")
  
  log_info("\n✅ Done.")
}

# Run if executed as script
if (sys.nframe() == 0) main(cfg)

