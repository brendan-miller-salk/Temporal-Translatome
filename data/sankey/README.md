# Sankey Diagram Data

This directory contains CSV files used to generate Sankey diagrams showing the classification pipeline for microproteins.

## File Naming Convention

Files are named according to classification path:
```
{SwissProt_status}_{DESeq_status}_{Enrichment}_{Localization}.csv
```

### Categories

**Swiss-Prot Status:**
- `Swiss-Prot`: Matches canonical Swiss-Prot protein
- `Non-Swiss-Prot`: Novel/non-canonical sequence

**DESeq Status:**
- `DESeqTrue`: Significantly differentially expressed (padj < 0.05)
- `DESeqFalse`: Not significantly differentially expressed

**Tissue Enrichment:**
- `Muscle-Enriched`: Tau score indicates muscle specificity
- `Non-Muscle-Enriched`: Broadly expressed or other tissue-specific

**Predicted Localization:**
- `Mitoch`: Mitochondrial localization predicted
- `Secreted`: Secreted/extracellular predicted
- `Other`: Cytoplasmic or other localization

## Example Files

| File | Description |
|------|-------------|
| `Non-Swiss-Prot_DESeqTrue_Muscle-Enriched_Mitoch.csv` | Novel, exercise-responsive, muscle-specific, mitochondrial microproteins |
| `Swiss-Prot_DESeqFalse_Non-Muscle-Enriched_Other.csv` | Canonical, stable expression, broadly expressed, cytoplasmic |

## Usage

These files are inputs for `code/figures/figure4_sankey.R` to generate the classification Sankey diagram. These may be used to prioritize functional experiments.
