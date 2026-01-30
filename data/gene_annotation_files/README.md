# Gene Annotation Files

GTF annotation files containing microprotein and smORF coordinates for use by the research community.

These files can be used directly with standard bioinformatics tools (featureCounts, bedtools, IGV, etc.) to quantify or visualize microprotein expression in your own datasets.

## Files

| File | Description |
|------|-------------|
| `Swiss-Prot.gtf` | Canonical Swiss-Prot annotated genes |
| `ShortStop_predictions.gtf` | ShortStop ML-predicted microproteins |
| `RiboSeq_evidence.gtf` | Transcripts with Ribo-seq translation evidence |
| `Mitochondria_predicted.gtf` | Predicted mitochondrial-localized microproteins |
| `Secreted_predicted.gtf` | Predicted secreted microproteins |
| `SEER_proteogenomics.gtf` | Mass spec-validated microproteins |
| `Exercise_differential_expression.csv` | Exercise-responsive genes for circos highlighting |

## GTF Format

Standard GTF format with additional attributes:
- `gene_id`: Unique identifier
- `transcript_id`: Transcript identifier  
- `gene_name`: Symbol (if available)
- `gene_biotype`: Classification (smORF type)

## Usage

These files can be used for omics analysis pipelines and genome visualization tools.
