#!/bin/bash
#===============================================================================
# STAR: RNA-seq Alignment for PANDA Mouse Atlas
#===============================================================================
# This script documents the STAR alignment parameters used to align
# bulk RNA-seq reads to the mouse genome for the PANDA microprotein atlas.
#
# Tool: STAR (Spliced Transcripts Alignment to a Reference)
# Reference: Dobin A, et al. (2013) STAR: ultrafast universal RNA-seq aligner.
#            Bioinformatics 29(1):15-21. PMID: 23104886
#===============================================================================

#-------------------------------------------------------------------------------
# EXACT COMMAND USED
#-------------------------------------------------------------------------------

STAR \
    --readFilesIn sample_R1.fastq.gz sample_R2.fastq.gz \
    --genomeDir STAR_genome_index \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMstrandField intronMotif \
    --sjdbOverhang 99 \
    --runThreadN 12 \
    --outFileNamePrefix output_prefix \
    --sjdbScore 2 \
    --outFilterMultimapNmax 2 \
    --outFilterScoreMinOverLread 0.25 \
    --outFilterMatchNminOverLread 0.25 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --alignSJoverhangMin 12 \
    --alignSJDBoverhangMin 1 \
    --readFilesCommand zcat

#-------------------------------------------------------------------------------
# PARAMETER EXPLANATION
#-------------------------------------------------------------------------------

# INPUT/OUTPUT
# --readFilesIn             Paired-end FASTQ files (R1 and R2)
# --genomeDir               Path to STAR genome index directory
# --outSAMtype              Output sorted BAM file
# --outSAMstrandField       Include strand information for downstream analysis
# --outFileNamePrefix       Prefix for output files
# --readFilesCommand zcat   Decompress gzipped FASTQ on-the-fly

# SPLICE JUNCTION DATABASE
# --sjdbOverhang            Read length - 1 (for 100bp reads, use 99)
# --sjdbScore               Score for annotated splice junctions

# MULTIMAPPING
# --outFilterMultimapNmax   Maximum number of loci a read can map to (2)
#                           Allows some multi-mapping for repetitive regions

# READ FILTERING
# --outFilterScoreMinOverLread     Minimum alignment score as fraction of read length
# --outFilterMatchNminOverLread    Minimum matched bases as fraction of read length
# --outFilterMismatchNoverReadLmax Maximum mismatches as fraction of read length

# ALIGNMENT PARAMETERS
# --alignIntronMin          Minimum intron size (20 bp)
# --alignIntronMax          Maximum intron size (1,000,000 bp)
# --alignMatesGapMax        Maximum gap between paired-end mates
# --alignSJoverhangMin      Minimum overhang for unannotated splice junctions
# --alignSJDBoverhangMin    Minimum overhang for annotated splice junctions

# THREADING
# --runThreadN              Number of threads for parallel processing

#-------------------------------------------------------------------------------
# INPUT FILES
#-------------------------------------------------------------------------------
# FASTQ files:    Paired-end Illumina RNA-seq reads (gzipped)
# Genome index:   STAR index built from GRCm39 primary assembly
#                 with GENCODE vM36 gene annotations
#                 Built with --sjdbOverhang matching read length - 1

#-------------------------------------------------------------------------------
# OUTPUT FILES
#-------------------------------------------------------------------------------
# *Aligned.sortedByCoord.out.bam   Coordinate-sorted BAM alignment
# *Log.final.out                   Alignment statistics summary
# *SJ.out.tab                      Detected splice junctions
# *Log.out                         Detailed run log

#-------------------------------------------------------------------------------
# GENOME INDEX GENERATION (for reference)
#-------------------------------------------------------------------------------
# STAR --runMode genomeGenerate \
#     --genomeDir STAR_genome_index \
#     --genomeFastaFiles GRCm39.primary_assembly.genome.fa \
#     --sjdbGTFfile gencode.vM36.annotation.gtf \
#     --sjdbOverhang 99 \
#     --runThreadN 12
