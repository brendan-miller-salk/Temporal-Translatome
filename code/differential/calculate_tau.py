#!/usr/bin/env python3
"""
calculate_tau.py - Tissue Specificity (Tau) Score Calculation

Calculates the tau tissue-specificity index for gene expression data.
Tau ranges from 0 (ubiquitous) to 1 (tissue-specific).

PANDA Atlas Methods:
Miller et al. (2025) - PANDA: Pan-tissue Atlas of Novel Differential microprotein Analysis

===============================================================================
BACKGROUND
===============================================================================
The tau index quantifies tissue specificity of gene expression. It was first
described by Yanai et al. (2005) "Genome-wide midrange transcription profiles 
reveal expression level relationships in human tissue specification."

Formula:
                      n
                     Σ (1 - x̂ᵢ)
                    i=1
            τ = ─────────────────
                     n - 1

Where:
- n = number of tissues
- x̂ᵢ = expression in tissue i normalized to max expression (xᵢ / max(x))

Interpretation:
- τ = 0: Gene is expressed equally across all tissues (housekeeping)
- τ = 1: Gene is expressed in only one tissue (tissue-specific)
- τ > 0.8: Generally considered tissue-specific (Kryuchkova-Mostacci & Robinson-Rechavi, 2017)
- τ > 0.9: Highly tissue-specific

===============================================================================
WORKFLOW
===============================================================================
1. Load raw count matrix from featureCounts/HTSeq
2. Load sample metadata with tissue annotations
3. Collapse tissue replicates (e.g., muscle subtypes → "Muscle")
4. Calculate CPM (Counts Per Million) normalization
5. Compute tau score for each gene
6. Identify tau_driver (tissue with highest expression)
7. Generate log2(CPM+1) matrix for downstream analysis

===============================================================================
PREREQUISITES
===============================================================================
- Python 3.7+
- pandas
- numpy

Input files:
- counts.csv: Raw count matrix (genes × samples)
- metadata.csv: Sample metadata with tissue annotations

===============================================================================
USAGE
===============================================================================
# As a module
from calculate_tau import calculate_tau_scores
tau_df = calculate_tau_scores(counts_file, metadata_file, output_file)

# As a script
python calculate_tau.py \\
    --counts combined_counts.csv \\
    --metadata sample_table.csv \\
    --output tau_scores.csv \\
    --method RNASeq

===============================================================================
OUTPUT
===============================================================================
CSV file with columns:
- gene_id: Gene/transcript identifier
- Tau: Tau score (0-1)
- tau_driver: Tissue with highest expression
- log2CPM_<tissue>: Log2(CPM+1) for each tissue

===============================================================================
"""

import pandas as pd
import numpy as np
import argparse
import sys


# =============================================================================
# TISSUE COLLAPSE MAPPING
# =============================================================================

# Define tissue groups for collapsing (customize as needed)
MUSCLE_TISSUES = [
    "Biceps femoris", 
    "Diaphragm", 
    "EDL", 
    "Gastrocnemius", 
    "Quadriceps",
    "Semimembranosus", 
    "Semitendinosus", 
    "Soleus", 
    "Tibialis anterior",
    "EDLTendon", 
    "TATendon", 
    "GasTendon"
]


def collapse_tissue(tissue, muscle_tissues=MUSCLE_TISSUES):
    """
    Collapse tissue subtypes into broader categories.
    
    Parameters:
    -----------
    tissue : str
        Original tissue name
    muscle_tissues : list
        List of tissue names to collapse into "Muscle"
    
    Returns:
    --------
    str : Collapsed tissue name
    """
    if tissue in muscle_tissues:
        return "Muscle"
    return tissue


# =============================================================================
# TAU SCORE CALCULATION
# =============================================================================

def calculate_tau(expression_row):
    """
    Calculate tau tissue-specificity index for a single gene.
    
    The tau index ranges from 0 (ubiquitous expression) to 1 (tissue-specific).
    
    Formula:
        τ = Σ(1 - x̂ᵢ) / (n - 1)
        where x̂ᵢ = xᵢ / max(x)
    
    Parameters:
    -----------
    expression_row : pd.Series
        Expression values across tissues (CPM normalized)
    
    Returns:
    --------
    float : Tau score between 0 and 1, or NaN if max expression is 0
    
    References:
    -----------
    Yanai et al. (2005) Bioinformatics 21(5):650-659
    Kryuchkova-Mostacci & Robinson-Rechavi (2017) Brief Bioinform 18(2):205-214
    """
    max_expr = expression_row.max()
    
    # Return NaN for unexpressed genes
    if max_expr == 0 or np.isnan(max_expr):
        return np.nan
    
    # Normalize to maximum expression
    normalized = expression_row / max_expr
    
    # Calculate tau
    n = len(expression_row)
    tau = (1 - normalized).sum() / (n - 1)
    
    # Clamp to [0, 1] for numerical stability
    return min(1.0, max(0.0, tau))


def get_tau_driver(expression_row):
    """
    Identify the tissue driving expression (highest expression).
    
    Parameters:
    -----------
    expression_row : pd.Series
        Expression values across tissues
    
    Returns:
    --------
    str : Name of tissue with highest expression, or NaN
    """
    if expression_row.max() > 0:
        return expression_row.idxmax()
    return np.nan


# =============================================================================
# MAIN WORKFLOW
# =============================================================================

def load_and_preprocess_counts(counts_file, metadata_file, method_filter="RNASeq"):
    """
    Load count matrix and metadata, then align samples.
    
    Parameters:
    -----------
    counts_file : str
        Path to raw counts CSV (genes × samples)
    metadata_file : str
        Path to sample metadata CSV
    method_filter : str
        Filter metadata by Method column (e.g., "RNASeq")
    
    Returns:
    --------
    tuple : (counts DataFrame, metadata DataFrame)
    """
    print(f"Loading counts from {counts_file}...")
    counts = pd.read_csv(counts_file, index_col=0)
    
    # Clean column names (remove alignment software suffixes)
    counts.columns = counts.columns.str.replace(
        r"Aligned\.sortedByCoord\.out\.bam", "", regex=True
    )
    
    print(f"Loading metadata from {metadata_file}...")
    metadata = pd.read_csv(metadata_file)
    
    # Filter by method if specified
    if method_filter and "Method" in metadata.columns:
        metadata = metadata[metadata["Method"] == method_filter]
        print(f"Filtered to {len(metadata)} {method_filter} samples")
    
    # Clean filename column
    if "FileName" in metadata.columns:
        metadata["FileName"] = metadata["FileName"].str.replace(
            r"_R[12]_001\.fastq\.gz|_R1_001\.sorted\.uniqueAlign\.bam", "", regex=True
        )
        metadata = metadata.drop_duplicates(subset="FileName")
    
    # Align samples between counts and metadata
    if "FileName" in metadata.columns:
        common_samples = [col for col in counts.columns if col in metadata["FileName"].values]
        counts = counts[common_samples]
        metadata = metadata.set_index("FileName").loc[common_samples]
        print(f"Aligned {len(common_samples)} common samples")
    
    return counts, metadata


def collapse_tissues(counts, metadata, tissue_column="Tissue"):
    """
    Collapse biological replicates by tissue category.
    
    Parameters:
    -----------
    counts : pd.DataFrame
        Raw count matrix (genes × samples)
    metadata : pd.DataFrame
        Sample metadata with tissue annotations
    tissue_column : str
        Column name containing tissue information
    
    Returns:
    --------
    pd.DataFrame : Collapsed count matrix (genes × tissues)
    """
    # Apply tissue collapsing
    if "Tissue_Collapsed" not in metadata.columns:
        metadata["Tissue_Collapsed"] = metadata[tissue_column].apply(collapse_tissue)
    
    # Rename columns to tissue names
    counts.columns = metadata["Tissue_Collapsed"].values
    
    # Collapse by taking mean across replicates
    collapsed = counts.groupby(by=counts.columns, axis=1).mean()
    
    print(f"Collapsed to {len(collapsed.columns)} tissue categories:")
    for tissue in sorted(collapsed.columns):
        print(f"  - {tissue}")
    
    return collapsed


def calculate_cpm(counts):
    """
    Calculate Counts Per Million (CPM) normalization.
    
    CPM = (raw_count / library_size) × 10^6
    
    Parameters:
    -----------
    counts : pd.DataFrame
        Raw or collapsed count matrix
    
    Returns:
    --------
    pd.DataFrame : CPM-normalized matrix
    """
    library_sizes = counts.sum(axis=0)
    cpm = (counts / library_sizes) * 1e6
    return cpm


def calculate_tau_scores(counts_file, metadata_file, output_file=None, 
                         method_filter="RNASeq", tissue_column="Tissue"):
    """
    Main function to calculate tau scores from raw counts.
    
    Parameters:
    -----------
    counts_file : str
        Path to raw counts CSV
    metadata_file : str
        Path to sample metadata CSV
    output_file : str, optional
        Path to save output CSV
    method_filter : str
        Filter for Method column in metadata
    tissue_column : str
        Column name for tissue annotations
    
    Returns:
    --------
    pd.DataFrame : Tau scores with log2CPM values
    """
    # Load and preprocess
    counts, metadata = load_and_preprocess_counts(
        counts_file, metadata_file, method_filter
    )
    
    # Collapse tissues
    collapsed = collapse_tissues(counts, metadata, tissue_column)
    
    # Calculate CPM
    print("Calculating CPM normalization...")
    cpm = calculate_cpm(collapsed)
    
    # Calculate tau scores
    print("Computing tau scores...")
    tau_results = pd.DataFrame({
        "Tau": cpm.apply(calculate_tau, axis=1),
        "tau_driver": cpm.apply(get_tau_driver, axis=1)
    })
    
    # Add log2(CPM+1) values
    log2_cpm = np.log2(cpm + 1)
    log2_cpm.columns = [f"log2CPM_{col}" for col in log2_cpm.columns]
    
    # Merge results
    tau_results = tau_results.merge(log2_cpm, left_index=True, right_index=True)
    
    # Summary statistics
    print("\n" + "="*60)
    print("TAU SCORE SUMMARY")
    print("="*60)
    print(f"Total genes: {len(tau_results)}")
    print(f"Genes with τ > 0.8 (tissue-specific): {(tau_results['Tau'] > 0.8).sum()}")
    print(f"Genes with τ > 0.9 (highly specific): {(tau_results['Tau'] > 0.9).sum()}")
    print(f"Genes with τ < 0.2 (housekeeping): {(tau_results['Tau'] < 0.2).sum()}")
    print("\nTau driver distribution:")
    print(tau_results['tau_driver'].value_counts().head(10))
    
    # Save results
    if output_file:
        tau_results.to_csv(output_file)
        print(f"\nResults saved to {output_file}")
    
    return tau_results


# =============================================================================
# COMMAND LINE INTERFACE
# =============================================================================

if __name__ == "__main__":
    print("""
 ___________  _   _   _ 
    |   /\   | | | | | |
    |  /--\  |_| |_| |_|
                        
    Tissue Specificity Score Calculator
    PANDA Atlas Pipeline
    """)
    
    parser = argparse.ArgumentParser(
        description="Calculate tau tissue-specificity scores from RNA-seq data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python calculate_tau.py --counts counts.csv --metadata samples.csv --output tau.csv
  
  python calculate_tau.py \\
      --counts combined_counts.csv \\
      --metadata sample_table.csv \\
      --output tau_scores.csv \\
      --method RNASeq \\
      --tissue Tissue

Interpretation:
  τ = 0:   Ubiquitous (housekeeping gene)
  τ > 0.8: Tissue-specific
  τ = 1:   Expressed in single tissue only
        """
    )
    
    parser.add_argument(
        "--counts", required=True,
        help="Path to raw count matrix CSV (genes × samples)"
    )
    parser.add_argument(
        "--metadata", required=True,
        help="Path to sample metadata CSV"
    )
    parser.add_argument(
        "--output", required=True,
        help="Output path for tau scores CSV"
    )
    parser.add_argument(
        "--method", default="RNASeq",
        help="Filter metadata by Method column (default: RNASeq)"
    )
    parser.add_argument(
        "--tissue", default="Tissue",
        help="Column name for tissue annotations (default: Tissue)"
    )
    
    args = parser.parse_args()
    
    try:
        tau_df = calculate_tau_scores(
            counts_file=args.counts,
            metadata_file=args.metadata,
            output_file=args.output,
            method_filter=args.method,
            tissue_column=args.tissue
        )
        print("\nTau calculation complete!")
        sys.exit(0)
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


# =============================================================================
# EXAMPLE USAGE (for documentation)
# =============================================================================
"""
# Example 1: Command line usage
python calculate_tau.py \\
    --counts combined_counts.csv \\
    --metadata panda_sample_table.csv \\
    --output tau_scores.csv \\
    --method RNASeq

# Example 2: Python module usage
from calculate_tau import calculate_tau_scores, calculate_tau

# Full pipeline
tau_df = calculate_tau_scores(
    counts_file="combined_counts.csv",
    metadata_file="sample_table.csv",
    output_file="tau_scores.csv"
)

# Filter for tissue-specific genes
tissue_specific = tau_df[tau_df['Tau'] > 0.9]
print(f"Found {len(tissue_specific)} highly tissue-specific genes")

# Get muscle-enriched genes
muscle_enriched = tau_df[tau_df['tau_driver'] == 'Muscle']

# Example 3: Calculate tau for custom expression data
import pandas as pd
import numpy as np

expression = pd.Series({
    'Brain': 100,
    'Heart': 5,
    'Liver': 2,
    'Muscle': 1,
    'Kidney': 3
})

tau = calculate_tau(expression)
print(f"Tau = {tau:.3f}")  # High tau indicates tissue-specificity (brain-enriched)

# Example 4: Interpreting results
# τ ~ 0.95: Expressed almost exclusively in one tissue
# τ ~ 0.75: Enriched in one tissue but detected in others  
# τ ~ 0.50: Moderately variable expression
# τ ~ 0.10: Near-ubiquitous expression (housekeeping)
"""
