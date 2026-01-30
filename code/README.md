# Analysis Code

This directory contains scripts used for all analyses. 

## Directory Structure

```
code/
├── alignment/           # Read alignment and quantification
│   ├── star_alignment.sh
│   ├── featureCounts_rnaseq.sh
│   ├── featureCounts_riboseq.sh
│   └── cutadapt_trim.sh
├── differential/        # Differential expression analysis
│   └── quads_deseq_pipeline.R
└── figures/             # Figure generation
    ├── figure1_summary.R
    ├── figure2_tau.R
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

1. **Adapter Trimming** (`cutadapt_trim.sh`)
   - Remove 3' adapters and low-quality bases from Ribo-seq reads

2. **Alignment** (`star_alignment.sh`)
   - Align reads to mouse genome (GRCm39) with STAR

3. **Quantification** (`featureCounts_*.sh`)
   - Count reads on GENCODE annotations and novel smORFs
   - Two-pass strategy to avoid double-counting

4. **Differential Expression** (`quads_deseq_pipeline.R`)
   - DESeq2 analysis with Sex × Exercise × Time design
   - Contrast extraction for each comparison

5. **Visualization** (`figures/`)
   - Generate publication-ready figures
