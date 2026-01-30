# Supplementary Tables

This directory contains supplementary data tables for the PANDA atlas manuscript.

## Tables

| File | Description |
|------|-------------|
| `S1_figure1_summary.csv` | Summary statistics for smORF types and evidence integration |
| `S2_figure2_tau_output.csv` | Tau tissue-specificity scores for all microproteins |
| `S3_figure3_exercise_volcano.csv` | DESeq2 results for exercise (RW vs Sed) comparison |
| `S4_figure3_sex_volcano.csv` | DESeq2 results for sex (F vs M) comparison |
| `S5_figure3_time_volcano.csv` | DESeq2 results for time-of-day (ZT4 vs ZT16) comparison |
| `S6_figure4_sankey_paths.csv` | Classification paths for Sankey diagram |
| `S7_quad_rpf_swissprot.csv` | Ribo-seq DE results for Swiss-Prot/Swiss-Prot-MP annotated proteins, Pilot analysis, n = 2 per group (RW vs Sed) |
| `S8_quad_rpf_non_swissprot.csv` | Ribo-seq DE results for non-Swiss-Prot microproteins, Pilot analysis, n = 2 per group (RW vs Sed) |

## Column Descriptions

### S7_quad_rpf_swissprot.csv

| Column | Description |
|--------|-------------|
| `UCSC_link` | Clickable link to UCSC Genome Browser for the gene symbol |
| `sequence` | Amino acid sequence of the protein |
| `gene_id` | Unique gene/transcript identifier with Swiss-Prot annotation |
| `log2FC` | Log2 fold change (RW vs Sed, Ribo-seq) |
| `pvalue` | Raw p-value for differential expression |
| `padj` | Benjamini-Hochberg adjusted p-value (FDR) |

### S8_fquad_rpf_non_swissprot.csv

| Column | Description |
|--------|-------------|
| `UCSC_link` | Clickable link to UCSC Genome Browser for the locus |
| `sequence` | Amino acid sequence of the microprotein |
| `gene_id` | Unique gene/transcript identifier |
| `gene_name` | Gene symbol (if available) |
| `smORF_type` | Classification type of the small open reading frame |
| `log2FC` | Log2 fold change (RW vs Sed, Ribo-seq) |
| `pvalue` | Raw p-value for differential expression |
| `padj` | Benjamini-Hochberg adjusted p-value (FDR) |

## Column Descriptions (continued)

### DESeq2 Output Tables (S3-S5)

| Column | Description |
|--------|-------------|
| `gene_id` | Unique gene/transcript identifier |
| `gene_symbol` | Gene symbol (if available) |
| `baseMean` | Mean normalized counts across all samples |
| `log2FoldChange` | Log2 fold change between conditions |
| `lfcSE` | Standard error of log2 fold change |
| `pvalue` | Raw p-value |
| `padj` | Benjamini-Hochberg adjusted p-value |

### Tau Scores (S2)

| Column | Description |
|--------|-------------|
| `gene_id` | Unique identifier |
| `tau` | Tissue specificity score (0 = ubiquitous, 1 = tissue-specific) |
| `max_tissue` | Tissue with highest expression |
| `tissue_*` | Normalized expression in each tissue |
