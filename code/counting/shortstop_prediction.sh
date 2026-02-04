#!/bin/bash
#===============================================================================
# ShortStop: Machine Learning Microprotein Discovery for PANDA Mouse Atlas
#===============================================================================
# This script documents the ShortStop analysis used to predict novel
# microproteins from Ribo-seq data in the PANDA mouse muscle atlas.
#
# Tool: ShortStop
# Reference: Miller B, de Souza EV, Pai VJ, Kim H, Vaughan JM, Lau CJ,
#            Diedrich JK, Saghatelian A. (2025) ShortStop: a machine learning
#            framework for microprotein discovery. BMC Methods.
#            PMID: 40756675 | DOI: 10.1186/s44330-025-00037-4
#
#===============================================================================

#-------------------------------------------------------------------------------
# STEP 1: Feature extraction
#-------------------------------------------------------------------------------
# Extract sequence and structural features from candidate smORFs

python ShortStop.py feature_extract \
    --genome genome.fa \
    --putative_smorfs_gtf candidate_smorfs.gtf \
    --outdir features_output/

#-------------------------------------------------------------------------------
# FEATURE EXTRACTION PARAMETERS
#-------------------------------------------------------------------------------
# --genome              Reference genome FASTA (GRCm39 primary assembly)
# --putative_smorfs_gtf GTF file with candidate smORF coordinates
# --outdir              Output directory for extracted features

#-------------------------------------------------------------------------------
# STEP 2: Run ShortStop prediction
#-------------------------------------------------------------------------------

python ShortStop.py predict \
    --genome genome.fa \
    --putative_smorfs_gtf candidate_smorfs.gtf \
    --orfs_to_be_predicted features_output/features.csv

#-------------------------------------------------------------------------------
# PREDICTION PARAMETERS
#-------------------------------------------------------------------------------
# --genome                        Reference genome FASTA (GRCm39 primary assembly)
# --putative_smorfs_gtf           GTF file with candidate smORF coordinates
# --orfs_to_be_predicted          Features extracted in Step 1

#-------------------------------------------------------------------------------
# INPUT FILES
#-------------------------------------------------------------------------------
# genome.fa:              GRCm39 primary assembly genome
# candidate_smorfs.gtf:   Candidate smORF coordinates from ORF calling

#-------------------------------------------------------------------------------
# OUTPUT FORMAT
#-------------------------------------------------------------------------------
# ShortStop outputs predictions with:
#   - ORF identifier
#   - SAM prediction 
#   - Confidence score (0-1)

#-------------------------------------------------------------------------------
# FILTERING CRITERIA USED IN PANDA ATLAS
#-------------------------------------------------------------------------------
# High-confidence microprotein predictions were filtered by:
#   1. SAM prediction
#   2. Confidence score > 0.5
#   3. Minimum CPM > 5 in at least one tissue

#-------------------------------------------------------------------------------
# INTEGRATION WITH PANDA ATLAS
#-------------------------------------------------------------------------------
# ShortStop predictions were combined with:
#   - Direct Ribo-seq evidence (tissue-specific translation)
#   - SEER proteogenomics (mass spectrometry validation)
#   - Swiss-Prot annotations (known microproteins)
