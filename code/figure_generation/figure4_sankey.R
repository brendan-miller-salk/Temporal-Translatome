# Panda Microprotein Sankey Plot
# Inspired by filtering criteria for microprotein evidence and differential expression

library(dplyr)
library(stringr)
library(networkD3)
library(readr)
library(tidyr)
library(htmltools)
library(htmlwidgets)
library(scales)
library(colorspace)
library(purrr)

# Load the data
data_path <- "/Users/brendanmiller/Library/CloudStorage/Box-Box/panda/dictionaries/panda_hits_deseq_deeploc.csv"
df <- read_csv(data_path)

table(df$proteogenomics)
# Apply the filtering logic as specified
filtered_df <- df %>%
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
table(filtered_df$has_seer)

# Create hierarchical splits for Sankey diagram (explicit node structure)
sankey_df <- filtered_df %>%
  rowwise() %>%
  mutate(
    # Check for ANY _padj < 0.05
    deseq_sig = {
      padj_vals <- c_across(ends_with("_padj"))
      if (length(padj_vals) == 0 || all(is.na(padj_vals))) {
        "False"
      } else {
        ifelse(any(padj_vals < 0.05, na.rm = TRUE), "True", "False")
      }
    },
    
    # Localization
    localization = case_when(
      str_detect(coalesce(Localizations, ""), "Mito") ~ "Mitoch",
      str_detect(coalesce(gene_ontology, ""), "mito") ~ "Mitoch",
      !is.na(Extracellular) & Extracellular > 0.5 ~ "Secreted",
      !is.na(`SP(Sec/SPI)`) & `SP(Sec/SPI)` > 0.5 ~ "Secreted",
      str_detect(coalesce(gene_ontology, ""), "extra|secrete|hormon") ~ "Secreted",
      TRUE ~ "Other"
    )
  ) %>%
  ungroup() %>%
  mutate(
    # Create explicit node hierarchy with arrows
    Root = "All Microproteins",
    Node1 = if_else(primary_smorf_type == "Swiss-Prot-MP", "Swiss-Prot", "Non-Swiss-Prot"),
    Node2 = paste("DESeq", deseq_sig),
    # Ensure we never generate an NA muscle bucket; default to Non-Muscle-Enriched
    Node3 = case_when(
      tau_driver == "Muscle" ~ "Muscle-Enriched",
      TRUE ~ "Non-Muscle-Enriched"
    ),
    # Build arrow-style combinations for terminal nodes
    combo_full = paste(Node1, Node2, Node3, localization, sep = " → "),
    Node4 = paste(Node2, "→", localization)
  ) %>%
  filter(!is.na(localization))

# Debug: Check filtering impact
cat("\nFiltering summary:\n")
cat("After initial filtering:", nrow(filtered_df), "\n")
cat("After localization assignment:", nrow(sankey_df), "\n")

# Build links between consecutive node levels
links_raw <- bind_rows(
  sankey_df %>% count(Root, Node1, name = "value") %>% rename(source = Root, target = Node1),
  sankey_df %>% count(Node1, Node2, name = "value") %>% rename(source = Node1, target = Node2),
  sankey_df %>% count(Node2, Node3, name = "value") %>% rename(source = Node2, target = Node3),
  sankey_df %>% count(Node3, Node4, name = "value") %>% rename(source = Node3, target = Node4)
)

# Build node counts
node_counts <- bind_rows(
  sankey_df %>% count(Root, name = "total") %>% rename(name = Root),
  sankey_df %>% count(Node1, name = "total") %>% rename(name = Node1),
  sankey_df %>% count(Node2, name = "total") %>% rename(name = Node2),
  sankey_df %>% count(Node3, name = "total") %>% rename(name = Node3),
  sankey_df %>% count(Node4, name = "total") %>% rename(name = Node4)
)

nodes <- data.frame(name = unique(c(links_raw$source, links_raw$target))) %>%
  left_join(node_counts, by = "name") %>%
  mutate(
    total = replace_na(total, 0),
    label = paste0(name, " (n=", total, ")")
  )

# Add grouping for color coding
# NOTE: Order matters! Check "Non-Swiss-Prot" before "Swiss-Prot" to avoid substring match
nodes$group <- case_when(
  grepl("DESeq True", nodes$name) ~ "DESeq_True",
  grepl("DESeq False", nodes$name) ~ "DESeq_False",
  grepl("Non-Swiss-Prot", nodes$name) ~ "NonSwiss",
  grepl("Swiss-Prot", nodes$name) ~ "Swiss",
  grepl("Non-Muscle-Enriched", nodes$name) ~ "NonMuscle",
  grepl("Muscle-Enriched", nodes$name) ~ "Muscle",
  grepl("^All Microproteins", nodes$name) ~ "Root",
  grepl("Mitoch", nodes$name) ~ "Mitoch",
  grepl("Secreted", nodes$name) ~ "Secreted",
  grepl("Other", nodes$name) ~ "Other_Loc",
  TRUE ~ "Other"
)

# Nature-quality color palette (colorblind-friendly, print-safe)
# Using volcano plot color theme for consistency
nodes$color <- case_when(
  nodes$group == "Root" ~ "#2D3142",
  nodes$group == "Swiss" ~ "#74a2b7",
  nodes$group == "NonSwiss" ~ "#FDE28D",
  nodes$group == "DESeq_True" ~ "#C44536",
  nodes$group == "DESeq_False" ~ "#626262",
  nodes$group == "Muscle" ~ "#b5b5b5",
  nodes$group == "NonMuscle" ~ "#313131",
  nodes$group == "Mitoch" ~ "#E8871E",
  nodes$group == "Secreted" ~ "#8B5A83",
  nodes$group == "Other_Loc" ~ "#B4C4AE",
  TRUE ~ "#fafafa"
)

# Build JS color scale - ensure group-color pairs are properly matched
group_color_pairs <- nodes %>%
  select(group, color) %>%
  distinct()

colour_map <- paste0(
  "d3.scaleOrdinal().domain([",
  paste0('"', group_color_pairs$group, '"', collapse = ", "),
  "]).range([",
  paste0('"', group_color_pairs$color, '"', collapse = ", "),
  "])"
)

links <- links_raw %>%
  mutate(
    IDsource = match(source, nodes$name) - 1,
    IDtarget = match(target, nodes$name) - 1
  )

# Publication-quality CSS styling for Nature
custom_css <- '
  <style>
    @import url("https://fonts.googleapis.com/css2?family=Source+Sans+Pro:wght@400;600&display=swap");
    
    /* Main text styling - clean sans-serif */
    text {
      font-family: "Source Sans Pro", "Helvetica Neue", Arial, sans-serif !important;
      fill: #1a1a1a !important;
      font-weight: 500 !important;
      font-size: 16px !important;
    }
    
    /* Node rectangles - subtle borders */
    .node rect {
      stroke: rgba(0,0,0,0.3) !important;
      stroke-width: 0.5px !important;
      rx: 2px;
      ry: 2px;
    }
    
    /* Link styling - semi-transparent with gradient feel */
    .link {
      stroke-opacity: 0.35 !important;
      transition: stroke-opacity 0.2s ease;
    }
    
    .link:hover {
      stroke-opacity: 0.7 !important;
    }
    
    /* Container styling */
    svg {
      background: #ffffff !important;
    }
  </style>
'

# Create the Sankey plot with publication-quality settings
sankey_plot <- sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = "IDsource",
  Target = "IDtarget",
  Value = "value",
  NodeID = "label",
  NodeGroup = "group",
  sinksRight = FALSE,
  fontSize = 18,
  fontFamily = "Source Sans Pro, Helvetica Neue, Arial, sans-serif",
  nodeWidth = 15,
  nodePadding = 25,
  margin = list(top = 30, right = 200, bottom = 30, left = 30),
  iterations = 100,
  colourScale = JS(colour_map),
  LinkGroup = NULL
)

# Prepend custom CSS
sankey_plot <- prependContent(sankey_plot, HTML(custom_css))

# Debug: Check the structure
cat("Links structure:\n")
print(links)
cat("\nNodes structure:\n")
print(nodes)

# Display the plot
sankey_plot

# Print summary statistics
cat("\nSummary:\n")
cat("Total proteins:", nrow(sankey_df), "\n")
cat("\nNode counts:\n")
print(node_counts)

# Save detailed combo counts (full path) so information is retained outside the simplified Sankey
write_csv(
  sankey_df %>% count(combo_full, name = "n") %>% arrange(desc(n)),
  "/Users/brendanmiller/Library/CloudStorage/Box-Box/panda/deseq/panda_microprotein_sankey_combos.csv"
)
cat("\nCombo details saved as 'panda_microprotein_sankey_combos.csv'\n")

# Save the plot as HTML with self-contained resources
htmlwidgets::saveWidget(
  sankey_plot, 
  "/Users/brendanmiller/Library/CloudStorage/Box-Box/panda/deseq/panda_microprotein_sankey.html",
  selfcontained = TRUE
)
cat("\nSankey plot saved as 'panda_microprotein_sankey.html'\n")

# Save as SVG for vector editing in Illustrator
fig_dir <- "/Users/brendanmiller/Library/CloudStorage/Box-Box/brendan_alan_shared_folder/Panda_paper/figures/figure4/panels/"
dir.create(fig_dir, showWarnings = FALSE)

# Export as SVG (vector format - editable in Illustrator)
if (requireNamespace("webshot2", quietly = TRUE)) {
  library(chromote)
  
  # Create a chromote session
  b <- ChromoteSession$new()
  
  # Navigate to the HTML file
  b$Page$navigate(paste0("file://", "/Users/brendanmiller/Library/CloudStorage/Box-Box/panda/deseq/panda_microprotein_sankey.html"))
  Sys.sleep(3)  # Wait for page to fully load
  
  # Extract the SVG element with proper namespace
  svg_content <- b$Runtime$evaluate("document.querySelector('svg').outerHTML")$result$value
  
  # Add proper XML/SVG headers for Illustrator compatibility
  svg_header <- '<?xml version="1.0" encoding="UTF-8"?>\n'
  
  # Add xmlns attribute if missing
  if (!grepl('xmlns=', svg_content)) {
    svg_content <- sub('<svg', '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"', svg_content)
  }
  
  # Combine header and content
  full_svg <- paste0(svg_header, svg_content)
  
  # Save as SVG file
  svg_path <- file.path(fig_dir, "sankey_microprotein_atlas.svg")
  writeLines(full_svg, svg_path, useBytes = TRUE)
  cat("SVG saved for Illustrator editing:", svg_path, "\n")
  
  # Close the browser session
  b$close()
  
  # Also save PDF (vector) for backup - Illustrator can also open this
  webshot2::webshot(
    "/Users/brendanmiller/Library/CloudStorage/Box-Box/panda/deseq/panda_microprotein_sankey.html",
    file.path(fig_dir, "sankey_microprotein_atlas.pdf"),
    vwidth = 1000,
    vheight = 600,
    delay = 1
  )
  cat("PDF saved (also vector-editable in Illustrator)\n")
} else {
  cat("Install webshot2 and chromote for SVG export: install.packages(c('webshot2', 'chromote'))\n")
}

cat("\n💡 TIP: If SVG has issues, open the HTML file in Chrome, right-click the diagram,\n")
cat("   select 'Inspect', find the <svg> element, right-click > Copy > Copy outerHTML,\n")
cat("   paste into a text file and save as .svg\n")

# === Export CSVs by full path ===
output_dir <- "/Users/brendanmiller/Library/CloudStorage/Box-Box/panda/deseq/sankey_files/"
dir.create(output_dir, showWarnings = FALSE)

sankey_export <- sankey_df %>%
  mutate(Path = combo_full) %>%
  filter(!is.na(Path)) %>%
  mutate(
    Path = gsub("→", "_", Path),
    Path = gsub("\\s+", "", Path)
  )

split(sankey_export, sankey_export$Path) %>%
  iwalk(~ {
    file_path <- file.path(output_dir, paste0(.y, ".csv"))
    write_csv(.x, file = file_path)
  })

cat("\nPath-specific CSVs exported to:", output_dir, "\n")

# === Export S6 Supplementary: All Sankey Paths in Single CSV ===
s6_supplementary <- sankey_df %>%
  dplyr::select(
    UCSC_link = CLICK_UCSC,
    sequence,
    gene_id = gene_ids,
    gene_name = gene_body,
    smORF_type = primary_smorf_type,
    pathway = combo_full,
    swiss_prot = Node1,
    deseq_significant = deseq_sig,
    muscle_enriched = Node3,
    localization
  )

write_csv(
  s6_supplementary,
  "/Users/brendanmiller/Library/CloudStorage/Box-Box/brendan_alan_shared_folder/Panda_paper/supplemental/S6_figure4_sankey_paths.csv"
)
cat("\nS6 Supplementary saved: S6_figure4_sankey_paths.csv with", nrow(s6_supplementary), "microproteins\n")

# === Volcano Plot: DESeq True + Mitochondria ===
library(ggplot2)

# Filter for DESeq True + Mitochondria
mito_deseq_df <- sankey_df %>%
  filter(deseq_sig == "True", localization == "Mitoch")

cat("\nMitochondrial DESeq True subset:", nrow(mito_deseq_df), "proteins\n")

# Prepare plot data - using Exercise_main as example (adjust as needed)
padj_thr <- 0.05
df_volcano <- mito_deseq_df %>%
  mutate(
    log2FC = Exercise_main_log2FoldChange,
    padj = Exercise_main_padj,
    status = if_else(primary_smorf_type == "Swiss-Prot-MP", "Swiss-Prot", "Noncanonical"),
    sig = padj < padj_thr & !is.na(padj)
  ) %>%
  filter(!is.na(log2FC) & !is.na(padj)) %>%
  # Arrange so Noncanonical is plotted first (behind Swiss-Prot)
  arrange(desc(status == "Noncanonical"))

# Custom theme
theme_musk <- function() {
  theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right"
    )
}

data_abs_max <- max(abs(df_volcano$log2FC), na.rm = TRUE)
padj_thr  <- 0.05                # Q-value threshold
x_cap_opt <- 0.30   
x_halfspan <- if (is.na(x_cap_opt)) data_abs_max else min(x_cap_opt, data_abs_max)
x_halfspan <- max(0.5, x_halfspan)  # never tighter than ±0.5
y_max <- max(-log10(df_volcano$padj), na.rm = TRUE)


# Create volcano plot
mito_volcano <- ggplot(df_volcano, aes(x = log2FC, y = -log10(padj))) +
  ggplot2::annotation_raster(
    raster = matrix(colorRampPalette(c("#cfe8ff", "white", "#f7b2b2"))(200), nrow = 1),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) +
  
  geom_point(
    aes(color = status, fill = status, alpha = sig),
    shape = 21, size = 3.8, stroke = 0.6
  ) +
  
  geom_hline(yintercept = -log10(padj_thr), linetype = "solid",
             color = "firebrick", linewidth = 1) +
  scale_fill_manual(values  = c("Swiss-Prot" = "#74a2b7", "Noncanonical" = "#FDE28D")) +
  scale_color_manual(values = c("Swiss-Prot" = "#5d91a9", "Noncanonical" = "#d3bb70")) +
  scale_alpha_manual(values = c(`TRUE` = 0.90, `FALSE` = 0.12), guide = "none") +
  labs(
    title = "Mitochondrial Microproteins",
    x = "Log2 Fold Change\nMain Effect Exercise",
    y = expression(-log[10]("Q-value"))
  ) +
  theme_minimal(base_family = "Roboto") +
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

# Print and save
print(mito_volcano)
ggsave(
  file.path("/Users/brendanmiller/Library/CloudStorage/Box-Box/brendan_alan_shared_folder/Panda_paper/figures/figure4/panels/mito_deseq_volcano.pdf"),
  mito_volcano,
  width = 4.75,
  height = 4.5,
  device = cairo_pdf
)
cat("\nMitochondrial volcano plot saved to figures/mito_deseq_volcano.pdf\n")

