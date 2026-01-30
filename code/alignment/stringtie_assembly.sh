#!/bin/bash
#===============================================================================
# StringTie: Transcript Assembly for PANDA Mouse Atlas
#===============================================================================
# This script documents the StringTie parameters used to assemble transcripts
# from aligned RNA-seq reads for the PANDA microprotein atlas.
#
# Tool: StringTie v2.1.6
# Reference: Pertea M, et al. (2015) StringTie enables improved reconstruction
#            of a transcriptome from RNA-seq reads. Nature Biotechnology 33:290-295.
#            PMID: 25690850
#===============================================================================

#-------------------------------------------------------------------------------
# EXACT COMMAND USED
#-------------------------------------------------------------------------------

stringtie \
    sample_Aligned.sortedByCoord.out.bam \
    -o sample.stringtie.gtf \
    -p 8 \
    -G gencode.vM33.annotation.gtf \
    --rf \
    -x chrM \
    -u \
    -v

#-------------------------------------------------------------------------------
# PARAMETER EXPLANATION
#-------------------------------------------------------------------------------

# INPUT/OUTPUT
# input.bam           Coordinate-sorted BAM file from STAR alignment
# -o <file>           Output GTF file with assembled transcripts
# -G <file>           Reference annotation GTF (guides assembly)

# LIBRARY TYPE
# --rf                Reverse-forward strand orientation (dUTP/Illumina TruSeq)
#                     First read maps to reverse strand of transcript

# FILTERING
# -x chrM             Exclude mitochondrial chromosome from assembly
#                     (mitochondrial transcripts can skew results)

# OUTPUT OPTIONS
# -u                  Disable multi-mapping correction
#                     (useful when multi-mappers already filtered)
# -v                  Verbose mode (print processing details)

# PERFORMANCE
# -p 8                Number of threads for parallel processing

#-------------------------------------------------------------------------------
# INPUT FILES
#-------------------------------------------------------------------------------
# BAM file:           STAR-aligned, coordinate-sorted RNA-seq reads
#                     (*Aligned.sortedByCoord.out.bam)
# Reference GTF:      GENCODE vM33 annotation for guided assembly
