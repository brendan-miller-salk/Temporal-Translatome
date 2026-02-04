# Analysis Code

This directory contains scripts used for all analyses. 

## Directory Structure

```
code/
├── alignment/              # Read alignment and assembly
│   ├── cutadapt_trim.sh
│   ├── star_alignment.sh
│   └── stringtie_assembly.sh
├── counting/               # Read quantification
│   ├── featureCounts_riboseq.sh
│   ├── featureCounts_rnaseq.sh
│   ├── make_unique_gtf_to_count_riboseq.sh
│   └── shortstop_prediction.sh
├── differential_expression/ # Differential expression analysis
│   ├── calculate_tau.py
│   └── quads_deseq_pipeline.R
└── figure_generation/      # Figure generation
    ├── figure1_summary.R
    ├── figure2_tau.r
    ├── figure3_deseq_volcano.R
    └── figure4_sankey.R
```

## Requirements

### R Packages
```r
# CRAN
install.packages(c("tidyverse", "ggplot2", "data.table", "ggrepel", "plotly"))

# Bioconductor
BiocManager::install(c("DESeq2", "biomaRt", "BiocParallel"))
```

### Command-Line Tools
- STAR v2.7.10+
- Subread (featureCounts) v2.0.2+
- cutadapt v3.5+
- StringTie v2.2+

## Pipeline Overview

1. **Adapter Trimming** (`alignment/cutadapt_trim.sh`)
   - Remove 3' adapters and low-quality bases from Ribo-seq reads

2. **Alignment** (`alignment/star_alignment.sh`)
   - Align reads to mouse genome (GRCm39) with STAR

3. **Assembly** (`alignment/stringtie_assembly.sh`)
   - Assemble transcripts with StringTie

4. **Quantification** (`counting/`)
   - `featureCounts_*.sh`: Count reads on GENCODE annotations and novel smORFs
   - `make_unique_gtf_to_count_riboseq.sh`: Prepare GTF for Ribo-seq counting
   - `shortstop_prediction.sh`: ShortStop ML-based smORF predictions

5. **Differential Expression** (`differential_expression/`)
   - `quads_deseq_pipeline.R`: DESeq2 analysis with Sex × Exercise × Time design
   - `calculate_tau.py`: Tissue specificity index calculation

6. **Visualization** (`figure_generation/`)
   - Generate publication-ready figures
