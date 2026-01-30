#!/bin/bash
#===============================================================================
# featureCounts: Bulk RNA-seq Read Quantification for PANDA Mouse Atlas
#===============================================================================
# This script documents the exact featureCounts command used to quantify
# gene/CDS expression from aligned BAM files in the PANDA microprotein atlas.
#
# Tool: featureCounts (Subread package v2.0.2)
# Reference: Liao Y, Smyth GK, Shi W. (2014) Bioinformatics 30(7):923-30
#===============================================================================

#-------------------------------------------------------------------------------
# EXACT COMMAND USED
#-------------------------------------------------------------------------------

featureCounts \
    -p \
    -O \
    -s 2 \
    -t CDS \
    -F GTF \
    -a annotation.gtf \
    -o sample_counts.txt \
    sample_Aligned.sortedByCoord.out.bam

#-------------------------------------------------------------------------------
# PARAMETER EXPLANATION
#-------------------------------------------------------------------------------
# -p              Paired-end mode: count fragments instead of individual reads
# -O              Multi-overlap: count reads overlapping multiple features
# -s 2            Reverse-stranded library (dUTP/Illumina TruSeq stranded)
# -t CDS          Count CDS features (coding sequences) for protein-level analysis
# -F GTF          Annotation format is GTF
# -a              GTF annotation file with gene/CDS coordinates
# -o              Output file for read counts

#-------------------------------------------------------------------------------
# INPUT FILES
#-------------------------------------------------------------------------------
# BAM files: STAR-aligned, coordinate-sorted (*_Aligned.sortedByCoord.out.bam)
# GTF file:  Custom annotation combining:
#            - GENCODE vM36 primary assembly (canonical genes)
#            - Ribo-seq validated microprotein ORFs
#            - ShortStop predicted smORFs
