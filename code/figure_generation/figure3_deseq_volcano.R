# === Volcano Plot (Clean, Muskified) ===
library(readr)
library(dplyr)
library(ggplot2)
library(stringr)


# ---- config ----
csv_path  <- "/Users/brendanmiller/Library/CloudStorage/Box-Box/panda/dictionaries/panda_hits_deseq.csv"
padj_thr  <- 0.05                # Q-value threshold
x_cap_opt <- 0.30                # optional symmetric cap; set to NA to auto from data

OUTDIR = '/Users/brendanmiller/Library/CloudStorage/Box-Box/panda/deseq/figures'

# ---- theme (Musk) ----
theme_musk <- function(base_family = "roboto") {
  theme_minimal(base_family = base_family) +
    theme(
      plot.title   = element_text(hjust = 0.5, size = 22),
      axis.title   = element_text(size = 18),
      axis.text    = element_text(size = 14, color = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(linetype = "dashed", linewidth = 0.3),
      panel.grid.major.y = element_line(linetype = "dashed", linewidth = 0.3),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.margin  = margin(8, 14, 8, 8),
      legend.position = "none",
      legend.title = element_blank()
    )
}

df_master <- read_csv(csv_path, show_col_types = FALSE) |>
  # ---- load & tidy ----
mutate(
  # unify primary_smorf_type
  primary_smorf_type = if_else(primary_smorf_type == "TrEMBL", "Noncanonical", primary_smorf_type),
  
  # guard against 0/NA padj for -log10 (keep numeric type)
  Exercise_main_padj = if_else(is.na(Exercise_main_padj) | Exercise_main_padj <= 0,
                               NA_real_, Exercise_main_padj),
  
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
 

# keep only relevant classes (drop others if present)
df_plot <- df_master |>
  rename(
    log2FC = Exercise_main_log2FoldChange,
    padj   = Exercise_main_padj
  ) |>
  mutate(
    sig = padj <= padj_thr,
    point_group = if_else(
      str_to_lower(primary_smorf_type) == "lncrna" &
        tau_driver == "Muscle" &
        str_detect(coalesce(shortstop_predictions, ""), "SAM"),
      "lncRNA",
      "other"
    )
  ) %>%
  mutate(point_group = factor(point_group, levels = c("other", "lncRNA"))) %>%
  arrange(point_group)

# ---- quick lookup: padj for a sequence in df_plot ----
sequence_of_interest <- "MRSGSCTSDWLPAVLSDCQTYSVLGGVGRDASSRGLPSDSHYPRLSISSFETEKRVLLRLSTSLDKEHSGHPGSIKYLAFLSLWNRLVYSVLFASERLIIKGSPTSSGTDRHPAPNPRQGLESSEDRDDAMREPFCGVAILC"
padj_lookup <- df_plot %>%
  filter(sequence == sequence_of_interest) %>%
  select(sequence, padj)
print(padj_lookup)

# ---- greatest absolute fold change (exercise) ----
max_fc_row <- df_plot %>%
  filter(padj < 0.001) %>%
  filter(!is.na(log2FC)) %>%
  filter(point_group == "lncRNA") %>%
  slice_max(order_by = abs(log2FC), n = 1, with_ties = FALSE) %>%
  select(sequence, log2FC, padj, primary_smorf_type)
print(max_fc_row)

# ---- symmetric x-range (optional cap, min ±0.5) ----
data_abs_max <- max(abs(df_plot$log2FC), na.rm = TRUE)
x_halfspan <- if (is.na(x_cap_opt)) data_abs_max else min(x_cap_opt, data_abs_max)
x_halfspan <- max(0.5, x_halfspan)  # never tighter than ±0.5
y_max <- max(-log10(df_plot$padj), na.rm = TRUE)

# ---- plot ----
exercise_volcano <-
  ggplot(df_plot, aes(x = log2FC, y = -log10(padj))) +
  ggplot2::annotation_raster(
    raster = matrix(colorRampPalette(c("#cfe8ff", "white", "#f7b2b2"))(200), nrow = 1),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) +
  geom_point(
    aes(alpha = point_group, size = point_group),
    color = "#847440", fill = "#FDE28D",
    shape = 21, stroke = 0.6
  ) +
  geom_point(
    data = df_plot %>%
      filter(sequence == "MRSGSCTSDWLPAVLSDCQTYSVLGGVGRDASSRGLPSDSHYPRLSISSFETEKRVLLRLSTSLDKEHSGHPGSIKYLAFLSLWNRLVYSVLFASERLIIKGSPTSSGTDRHPAPNPRQGLESSEDRDDAMREPFCGVAILC"),
    color = "black", fill = "#b22222",
    shape = 21, size = 10.6, alpha = 1, stroke = 0.8
  ) +
  geom_hline(yintercept = -log10(padj_thr), linetype = "solid",
             color = "firebrick", linewidth = 1) +
  scale_alpha_manual(values = c("lncRNA" = 1, "other" = 0.1), guide = "none") +
  scale_size_manual(values = c("lncRNA" = 7.0, "other" = 2.2), guide = "none") +
  labs(
    title = "Effect of Exercise",
    x = "Log2 Fold Change",
    y = expression(-log[10]("Q-value"))
  ) +
  #coord_cartesian(xlim = c(-x_halfspan, x_halfspan)) +
  theme_musk()

# print it
exercise_volcano
ggsave(file.path("/Users/brendanmiller/Library/CloudStorage/Box-Box/panda/deseq/figures/exercise_volcano.pdf"), exercise_volcano, width = 7, height = 3.5, device = cairo_pdf)

# === Save Exercise Volcano Supplementary ===
exercise_supplementary <- df_master %>%
  dplyr::select(
    UCSC_link = CLICK_UCSC,
    sequence,
    gene_id = gene_ids,
    gene_name = gene_body,
    smORF_type = primary_smorf_type,
    log2FC = Exercise_main_log2FoldChange,
    pvalue = Exercise_main_pvalue,
    padj = Exercise_main_padj
  )
write_csv(exercise_supplementary, "/Users/brendanmiller/Library/CloudStorage/Box-Box/brendan_alan_shared_folder/Panda_paper/supplemental/S3_figure3_exercise_volcano.csv")


# keep only relevant classes (drop others if present)

df_plot <- df_master |>
  rename(
    log2FC = Sex_Exercise_log2FoldChange,
    padj   = Sex_Exercise_padj
  ) |>
  mutate(
    sig = padj <= padj_thr,
    point_group = if_else(
      str_to_lower(primary_smorf_type) == "lncrna" &
        tau_driver == "Muscle" &
        str_detect(coalesce(shortstop_predictions, ""), "SAM"),
      "lncRNA",
      "other"
    )
  ) %>%
  mutate(point_group = factor(point_group, levels = c("other", "lncRNA"))) %>%
  arrange(point_group)

# ---- symmetric x-range (optional cap, min ±0.5) ----
data_abs_max <- max(abs(df_plot$log2FC), na.rm = TRUE)
x_halfspan <- if (is.na(x_cap_opt)) data_abs_max else min(x_cap_opt, data_abs_max)
x_halfspan <- max(0.5, x_halfspan)  # never tighter than ±0.5
y_max <- max(-log10(df_plot$padj), na.rm = TRUE)

# ---- plot ----
sex_volcano <-
  ggplot(df_plot, aes(x = log2FC, y = -log10(padj))) +
  ggplot2::annotation_raster(
    raster = matrix(colorRampPalette(c("#cfe8ff", "white", "#f7b2b2"))(200), nrow = 1),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) +
  geom_point(
    aes(alpha = point_group, size = point_group),
    color = "#847440", fill = "#FDE28D",
    shape = 21, stroke = 0.6
  ) +
  geom_point(
    data = df_plot %>%
      filter(sequence == "MVPVSFPLILVYGQTLVSHF"),
    color = "black", fill = "#b22222",
    shape = 21, size = 10.6, alpha = 1, stroke = 0.8
  ) +
  geom_hline(yintercept = -log10(padj_thr), linetype = "solid",
             color = "firebrick", linewidth = 1) +
  scale_alpha_manual(values = c("lncRNA" = 1, "other" = 0.1), guide = "none") +
  scale_size_manual(values = c("lncRNA" = 7.0, "other" = 2.2), guide = "none") +
  labs(
    title = "Effect of Biological Sex",
    x = "Log2 Fold Change",
    y = expression(-log[10]("Q-value"))
  ) +
  #coord_cartesian(xlim = c(-x_halfspan, x_halfspan)) +
  theme_musk()

# print it
sex_volcano
ggsave(file.path("/Users/brendanmiller/Library/CloudStorage/Box-Box/panda/deseq/figures/sex_volcano.pdf"), sex_volcano, width = 7, height = 3.5, device = cairo_pdf)

# === Save Sex Volcano Supplementary ===
sex_supplementary <- df_master %>%
  dplyr::select(
    UCSC_link = CLICK_UCSC,
    sequence,
    gene_id = gene_ids,
    gene_name = gene_body,
    smORF_type = primary_smorf_type,
    log2FC = Sex_Exercise_log2FoldChange,
    pvalue = Sex_Exercise_pvalue,
    padj = Sex_Exercise_padj
  )
write_csv(sex_supplementary, "/Users/brendanmiller/Library/CloudStorage/Box-Box/brendan_alan_shared_folder/Panda_paper/supplemental/S4_figure3_sex_volcano.csv")




# keep only relevant classes (drop others if present)
df_plot <- df_master |>
  rename(
    log2FC = Exercise_ZT_Interaction_log2FoldChange,
    padj   = Exercise_ZT_Interaction_padj
  ) |>
  mutate(
    sig = padj <= padj_thr,
    point_group = if_else(
      str_to_lower(primary_smorf_type) == "lncrna" &
        tau_driver == "Muscle" &
        str_detect(coalesce(shortstop_predictions, ""), "SAM"),
      "lncRNA",
      "other"
    )
  ) %>%
  mutate(point_group = factor(point_group, levels = c("other", "lncRNA"))) %>%
  arrange(point_group)

# ---- symmetric x-range (optional cap, min ±0.5) ----
data_abs_max <- max(abs(df_plot$log2FC), na.rm = TRUE)
x_halfspan <- if (is.na(x_cap_opt)) data_abs_max else min(x_cap_opt, data_abs_max)
x_halfspan <- max(0.5, x_halfspan)  # never tighter than ±0.5
y_max <- max(-log10(df_plot$padj), na.rm = TRUE)

# ---- plot ----
time_volcano <-
  ggplot(df_plot, aes(x = log2FC, y = -log10(padj))) +
  ggplot2::annotation_raster(
    raster = matrix(colorRampPalette(c("#cfe8ff", "white", "#f7b2b2"))(200), nrow = 1),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) +
  geom_point(
    aes(alpha = point_group, size = point_group),
    color = "#847440", fill = "#FDE28D",
    shape = 21, stroke = 0.6
  ) +
  geom_point(
    data = df_plot %>%
      filter(sequence == "MDSEASLGYQRRLCFKKKKKKRKLKWTNPQVLREWQGAYWALPTATMLPL"),
    color = "black", fill = "#b22222",
    shape = 21, size = 10.6, alpha = 1, stroke = 0.8
  ) +
  geom_hline(yintercept = -log10(padj_thr), linetype = "solid",
             color = "firebrick", linewidth = 1) +
  scale_alpha_manual(values = c("lncRNA" = 1, "other" = 0.1), guide = "none") +
  scale_size_manual(values = c("lncRNA" = 7.0, "other" = 2.2), guide = "none") +
  labs(
    title = "Effect of Time on Exercise",
    x = "Log2 Fold Change",
    y = expression(-log[10]("Q-value"))
  ) +
  #coord_cartesian(xlim = c(-x_halfspan, x_halfspan)) +
  theme_musk()

# print it
time_volcano
ggsave(file.path("/Users/brendanmiller/Library/CloudStorage/Box-Box/panda/deseq/figures/time_volcano.pdf"), time_volcano, width = 7, height = 3.5, device = cairo_pdf)

# === Save Time (Exercise x ZT Interaction) Volcano Supplementary ===
time_supplementary <- df_master %>%
  dplyr::select(
    UCSC_link = CLICK_UCSC,
    sequence,
    gene_id = gene_ids,
    gene_name = gene_body,
    smORF_type = primary_smorf_type,
    log2FC = Exercise_ZT_Interaction_log2FoldChange,
    pvalue = Exercise_ZT_Interaction_pvalue,
    padj = Exercise_ZT_Interaction_padj
  )
write_csv(time_supplementary, "/Users/brendanmiller/Library/CloudStorage/Box-Box/brendan_alan_shared_folder/Panda_paper/supplemental/S5_figure3_time_volcano.csv")




df_plot <- df_master |>
  rename(
    log2FC = Exercise_main_log2FoldChange,
    padj   = Exercise_main_padj
  ) |>
  mutate(
    sig = padj <= padj_thr,
    point_group = if_else(
      str_to_lower(primary_smorf_type) == "lncrna" &
        tau_driver == "Muscle" &
        str_detect(coalesce(shortstop_predictions, ""), "SAM"),
      "lncRNA",
      "other"
    )
  ) %>%
  mutate(point_group = factor(point_group, levels = c("other", "lncRNA"))) %>%
  arrange(point_group)

# ---- symmetric x-range (optional cap, min ±0.5) ----
data_abs_max <- max(abs(df_plot$log2FC), na.rm = TRUE)
x_halfspan <- if (is.na(x_cap_opt)) data_abs_max else min(x_cap_opt, data_abs_max)
x_halfspan <- max(0.5, x_halfspan)  # never tighter than ±0.5
y_max <- max(-log10(df_plot$padj), na.rm = TRUE)

# ---- plot ----
sex_volcano <-
  ggplot(df_plot, aes(x = log2FC, y = -log10(padj))) +
  ggplot2::annotation_raster(
    raster = matrix(colorRampPalette(c("#cfe8ff", "white", "#f7b2b2"))(200), nrow = 1),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) +
  geom_point(
    aes(alpha = point_group, size = point_group),
    color = "#847440", fill = "#FDE28D",
    shape = 21, stroke = 0.6
  ) +
  geom_point(
    data = df_plot %>%
      filter(sequence == "MVPVSFPLILVYGQTLVSHF"),
    color = "black", fill = "#b22222",
    shape = 21, size = 10.6, alpha = 1, stroke = 0.8
  ) +
  geom_hline(yintercept = -log10(padj_thr), linetype = "solid",
             color = "firebrick", linewidth = 1) +
  scale_alpha_manual(values = c("lncRNA" = 1, "other" = 0.1), guide = "none") +
  scale_size_manual(values = c("lncRNA" = 7.0, "other" = 2.2), guide = "none") +
  labs(
    title = "Effect of Biological Sex",
    x = "Log2 Fold Change",
    y = expression(-log[10]("Q-value"))
  ) +
  #coord_cartesian(xlim = c(-x_halfspan, x_halfspan)) +
  theme_musk()


# ============================================================
# === Combined Supplemental Figure: All Volcano Plots ===
# ============================================================
library(patchwork)

# Define all DESeq comparisons
deseq_comparisons <- list(
  list(name = "Exercise (Main)", log2fc = "Exercise_main_log2FoldChange", padj = "Exercise_main_padj"),
  list(name = "Time (Main)", log2fc = "Time_main_log2FoldChange", padj = "Time_main_padj"),
  list(name = "Sex (Main)", log2fc = "Sex_main_log2FoldChange", padj = "Sex_main_padj"),
  list(name = "Exercise × ZT", log2fc = "Exercise_ZT_Interaction_log2FoldChange", padj = "Exercise_ZT_Interaction_padj"),
  list(name = "Sex × Exercise", log2fc = "Sex_Exercise_log2FoldChange", padj = "Sex_Exercise_padj"),
  list(name = "Sex × Time", log2fc = "Sex_Time_log2FoldChange", padj = "Sex_Time_padj"),
  list(name = "Sex × Exercise × Time", log2fc = "Sex_Exercise_Time_log2FoldChange", padj = "Sex_Exercise_Time_padj"),
  list(name = "Exercise @ ZT4", log2fc = "Exercise_ZT4_log2FoldChange", padj = "Exercise_ZT4_padj"),
  list(name = "Exercise @ ZT16", log2fc = "Exercise_ZT16_log2FoldChange", padj = "Exercise_ZT16_padj"),
  list(name = "ZT4 Sedentary", log2fc = "ZT4_Sedentary_log2FoldChange", padj = "ZT4_Sedentary_padj"),
  list(name = "ZT4 Exercise", log2fc = "ZT4_Exercise_log2FoldChange", padj = "ZT4_Exercise_padj"),
  list(name = "Female @ ZT16", log2fc = "Female_ZT16_log2FoldChange", padj = "Female_ZT16_padj"),
  list(name = "Female @ ZT4", log2fc = "Female_ZT4_log2FoldChange", padj = "Female_ZT4_padj"),
  list(name = "Female Exercise", log2fc = "Female_Ex_log2FoldChange", padj = "Female_Ex_padj")
)

# Simplified theme for small multiples
theme_volcano_small <- function() {
  theme_minimal(base_family = "roboto") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 6, color = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linetype = "dashed", linewidth = 0.2),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.margin = margin(4, 4, 4, 4),
      legend.position = "none"
    )
}

# Function to create a single volcano plot
make_volcano <- function(df, log2fc_col, padj_col, title, padj_threshold = 0.05) {
  # Check if columns exist
  if (!log2fc_col %in% names(df) || !padj_col %in% names(df)) {
    return(ggplot() + theme_void() + labs(title = paste(title, "(N/A)")))
  }
  
  plot_df <- df %>%
    mutate(
      log2FC = .data[[log2fc_col]],
      padj = .data[[padj_col]]
    ) %>%
    filter(!is.na(log2FC) & !is.na(padj) & padj > 0)
  
  # Calculate symmetric x-axis limits centered at 0
  x_max <- max(abs(plot_df$log2FC), na.rm = TRUE)
  x_limit <- max(0.5, x_max)  # minimum ±0.5
  
  ggplot(plot_df, aes(x = log2FC, y = -log10(padj))) +
    ggplot2::annotation_raster(
      raster = matrix(colorRampPalette(c("#cfe8ff", "white", "#f7b2b2"))(200), nrow = 1),
      xmin = -x_limit, xmax = x_limit, ymin = -Inf, ymax = Inf
    ) +
    geom_point(
      color = "#847440", fill = "#FDE28D",
      shape = 21, stroke = 0.3, size = 1, alpha = 0.5
    ) +
    geom_hline(yintercept = -log10(padj_threshold), linetype = "solid",
               color = "firebrick", linewidth = 0.5) +
    coord_cartesian(xlim = c(-x_limit, x_limit)) +
    labs(title = title, x = "Log2FC", y = "-log10(Q)") +
    theme_volcano_small()
}

# Generate all volcano plots
volcano_plots <- lapply(deseq_comparisons, function(comp) {
  make_volcano(df_master, comp$log2fc, comp$padj, comp$name, padj_thr)
})

# Combine into a grid (4 rows × 4 cols = 16 panels, 14 used)
combined_volcano <- wrap_plots(volcano_plots, ncol = 4, nrow = 4)

# Save combined supplemental figure
ggsave(
  filename = "/Users/brendanmiller/Library/CloudStorage/Box-Box/brendan_alan_shared_folder/Panda_paper/supplemental/SFig4_volcano.pdf",
  plot = combined_volcano,
  width = 14, height = 12,
  device = cairo_pdf
)
cat("Saved combined volcano figure: SFig1_volcano.pdf\n")
