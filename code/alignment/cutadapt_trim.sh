#!/bin/bash
#===============================================================================
# cutadapt: Adapter Trimming for Ribo-seq FASTQ Files
#===============================================================================
# This script documents the exact cutadapt command used to remove adapter
# sequences and low-quality bases from ribosome profiling (Ribo-seq) reads.
#
# Tool: cutadapt v3.5
# Reference: Martin M. (2011) EMBnet.journal 17(1):10-12
#===============================================================================

#-------------------------------------------------------------------------------
# STEP 1: Trim adapter and filter low-quality reads
#-------------------------------------------------------------------------------
# Parameters:
#   -a ADAPTER    : 3' adapter sequence (identified via FastQC)
#   -m 15         : Discard reads shorter than 15 bp after trimming
#   -q 20         : Trim low-quality bases (Q < 20) from 3' end
#   -j 4          : Use 4 CPU threads

cutadapt \
    -a ADAPTER_SEQUENCE \
    -m 15 \
    -q 20 \
    -j 4 \
    -o sample_trimmed.fastq.gz \
    sample.fastq.gz
