#!/bin/bash
#===============================================================================
# featureCounts: Ribo-seq Read Quantification for PANDA Mouse Atlas
#===============================================================================
# This script documents the exact featureCounts commands used to quantify
# ribosome-protected fragment (RPF) counts for both canonical genes and
# novel microprotein ORFs in the PANDA atlas.
#
# Strategy: Two-pass counting to avoid double-counting overlapping features
#   1. Count reads on GENCODE canonical CDS
#   2. Count reads on unique (non-overlapping) smORF segments
#   3. Combine into unified count matrix
#
# Tool: featureCounts (Subread package v2.0.2)
# Reference: Liao Y, Smyth GK, Shi W. (2014) Bioinformatics 30(7):923-30
#===============================================================================

#-------------------------------------------------------------------------------
# STEP 1: Count GENCODE canonical CDS
#-------------------------------------------------------------------------------

featureCounts \
    -M \
    -O \
    -s 1 \
    -t CDS \
    -F GTF \
    -a GENCODEv36_cds.gtf \
    -o sample_gencode.txt \
    sample_sorted.uniqueAlign.bam

#-------------------------------------------------------------------------------
# STEP 2: Count unique non-overlapping smORF portions
#-------------------------------------------------------------------------------
# The unique GTF contains smORF CDS segments that do NOT overlap with any
# GENCODE annotation, ensuring no double-counting of shared genomic regions.

featureCounts \
    -M \
    -O \
    -s 1 \
    -t CDS \
    -F GTF \
    -a all_non_overlapping_portions.gtf \
    -o sample_unique.txt \
    sample_sorted.uniqueAlign.bam

#-------------------------------------------------------------------------------
# STEP 3: Combine count tables
#-------------------------------------------------------------------------------
# Merge GENCODE and unique smORF counts into a single matrix

# Extract counts (skip header, get gene_id and count columns)
tail -n +2 sample_gencode.txt | cut -f1,7 > sample_gencode_counts.txt
tail -n +2 sample_unique.txt | cut -f1,7 > sample_unique_counts.txt

# Concatenate into final output
echo -e "Geneid\tsample" > sample_combined.txt
cat sample_gencode_counts.txt sample_unique_counts.txt >> sample_combined.txt

#-------------------------------------------------------------------------------
# PARAMETER EXPLANATION
#-------------------------------------------------------------------------------
# -M              Count multi-mapping reads (important for repetitive smORFs)
# -O              Multi-overlap: count reads overlapping multiple features
# -s 1            Forward-stranded library (typical for Ribo-seq protocols)
# -t CDS          Count CDS features (coding sequences under translation)
# -F GTF          Annotation format is GTF
# -a              GTF annotation file

#-------------------------------------------------------------------------------
# INPUT FILES
#-------------------------------------------------------------------------------
# BAM files: Uniquely aligned, coordinate-sorted Ribo-seq reads
#            (*_sorted.uniqueAlign.bam)
#
# GTF files:
#   1. GENCODEv36_cds.gtf
#      - All CDS features from GENCODE vM36 primary assembly
#      - Canonical protein-coding gene annotations
#
#   2. all_non_overlapping_portions.gtf
#      - smORF CDS segments that do NOT overlap GENCODE annotations
#      - Derived by subtracting GENCODE coordinates from full smORF CDS
#      - Ensures reads are not double-counted between canonical and novel ORFs

#-------------------------------------------------------------------------------
# OUTPUT
#-------------------------------------------------------------------------------
# Combined count matrix with:
#   - GENCODE gene counts (canonical proteins)
#   - Unique smORF segment counts (novel microproteins)
#
# This two-GTF approach allows accurate quantification of novel smORFs
# while preserving canonical gene expression estimates.
