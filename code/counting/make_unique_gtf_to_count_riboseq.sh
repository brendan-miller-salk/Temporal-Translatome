#!/bin/bash
set -euo pipefail

# === Paths ===
REF_DIR="/Users/brendanmiller/Library/CloudStorage/Box-Box/panda/reference_databases"
GTF="$REF_DIR/rp3_riboseq_shortstop_panda_gencode.vM36.primary_assembly.basic.annotation.gtf"
OUTDIR="$REF_DIR/unique_gtf2fasta"
mkdir -p "$OUTDIR"

echo "🔍 Extracting non-overlapping portions of ALL GTF2FASTA CDS entries..."

# === Step 1: Extract CDS entries by source ===
echo "Step 1: Extracting CDS entries..."
awk -F'\t' '$2 == "GTF2FASTA" && $3 == "CDS"' "$GTF" > "$OUTDIR/gtf2fasta_cds.gtf"
awk -F'\t' '$2 != "GTF2FASTA" && $3 == "CDS"' "$GTF" > "$OUTDIR/GENCODEv36_cds.gtf"

echo "  GTF2FASTA CDS entries: $(wc -l < "$OUTDIR/gtf2fasta_cds.gtf")"
echo "  GENCODEv36 CDS entries: $(wc -l < "$OUTDIR/GENCODEv36_cds.gtf")"

# === Step 2: Convert to BED for bedtools ===
echo "Step 2: Converting to BED format..."
# Add line numbers to track original GTF entries
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $4-1, $5, NR, ".", $7}' "$OUTDIR/gtf2fasta_cds.gtf" > "$OUTDIR/gtf2fasta_cds.bed"
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $4-1, $5, ".", ".", $7}' "$OUTDIR/GENCODEv36_cds.gtf" > "$OUTDIR/GENCODEv36_cds.bed"

# === Step 2b: Extend GENCODE first/last CDS by 35nt per transcript ===
# This ensures we exclude GTF2FASTA regions within 35nt of GENCODE start/stop codons
# (where ribosome footprints from canonical genes might appear)
if [ -f "$OUTDIR/GENCODEv36_cds_flanked.bed" ]; then
    echo "Step 2b: Skipping (GENCODEv36_cds_flanked.bed already exists)..."
    echo "  GENCODE CDS (flanked): $(wc -l < "$OUTDIR/GENCODEv36_cds_flanked.bed")"
else
    echo "Step 2b: Extending GENCODE first/last CDS by 35nt..."

    awk -F'\t' 'BEGIN{OFS="\t"}
{
    # Extract transcript_id from attributes (BSD awk compatible)
    tid = $9
    gsub(/.*transcript_id "/, "", tid)
    gsub(/".*/, "", tid)
    
    # Store entry info
    n = ++count[tid]
    lines[tid, n] = $0
    chrs[tid, n] = $1
    starts[tid, n] = $4
    ends[tid, n] = $5
    strands[tid, n] = $7
}
END {
    for (tid in count) {
        n_entries = count[tid]
        strand = strands[tid, 1]
        
        # Find first CDS (smallest start on +, largest end on -)
        # Find last CDS (largest end on +, smallest start on -)
        first_idx = 1
        last_idx = 1
        
        for (i = 1; i <= n_entries; i++) {
            if (strand == "+") {
                if (starts[tid, i] < starts[tid, first_idx]) first_idx = i
                if (ends[tid, i] > ends[tid, last_idx]) last_idx = i
            } else {
                if (ends[tid, i] > ends[tid, first_idx]) first_idx = i
                if (starts[tid, i] < starts[tid, last_idx]) last_idx = i
            }
        }
        
        # Output BED with extensions for first/last CDS
        for (i = 1; i <= n_entries; i++) {
            bed_start = starts[tid, i] - 1  # Convert to 0-based
            bed_end = ends[tid, i]
            
            # Extend first CDS upstream by 35nt
            if (i == first_idx) {
                if (strand == "+") {
                    bed_start = (bed_start > 35) ? bed_start - 35 : 0
                } else {
                    bed_end = bed_end + 35
                }
            }
            
            # Extend last CDS downstream by 35nt
            if (i == last_idx) {
                if (strand == "+") {
                    bed_end = bed_end + 35
                } else {
                    bed_start = (bed_start > 35) ? bed_start - 35 : 0
                }
            }
            
            print chrs[tid, i], bed_start, bed_end, tid, ".", strand
        }
    }
}' "$OUTDIR/GENCODEv36_cds.gtf" > "$OUTDIR/GENCODEv36_cds_flanked.bed"

    echo "  GENCODE CDS (flanked): $(wc -l < "$OUTDIR/GENCODEv36_cds_flanked.bed")"
fi

# === Step 3: Extract non-overlapping portions (GTF2FASTA vs GENCODE+35nt+35nt) ===
echo "Step 3: Extracting non-overlapping CDS portions (excluding ±35nt around GENCODE start/stop)..."
bedtools subtract -s \
  -a "$OUTDIR/gtf2fasta_cds.bed" \
  -b "$OUTDIR/GENCODEv36_cds_flanked.bed" \
  > "$OUTDIR/gtf2fasta_nonoverlap_vs_gencode.bed"

echo "  Non-overlapping regions: $(wc -l < "$OUTDIR/gtf2fasta_nonoverlap_vs_gencode.bed")"
echo "  Original GTF2FASTA CDS: $(wc -l < "$OUTDIR/gtf2fasta_cds.bed")"

# === Step 4: Map non-overlapping BED segments back to GTF entries ===
echo "Step 4: Creating unique GTF entries for non-overlapping CDS segments..."

awk 'BEGIN{OFS="\t"}
NR==FNR {
    # Record BED info indexed by original line number (col 4)
    bed_chr[$4]=$1
    bed_start[$4]=$2+1
    bed_end[$4]=$3
    bed_strand[$4]=$6
    next
}
{
    line_num=FNR
    if (line_num in bed_chr) {
        print bed_chr[line_num], "GTF2FASTA", "CDS", \
              bed_start[line_num], bed_end[line_num], ".", \
              bed_strand[line_num], ".", $9
    }
}' "$OUTDIR/gtf2fasta_nonoverlap_vs_gencode.bed" "$OUTDIR/gtf2fasta_cds.gtf" \
> "$OUTDIR/gtf2fasta_unique_segments.gtf"

echo "  Created unique GTF entries: $(wc -l < "$OUTDIR/gtf2fasta_unique_segments.gtf")"

# === Step 5: Diagnostic — identify original CDS lines that split into multiple portions ===
echo "Step 5: Checking for split CDS entries..."
awk 'BEGIN{OFS="\t"}
{
    line_num=$4
    count[line_num]++
}
END {
    for (ln in count)
        if (count[ln]>1)
            print "Original CDS line", ln, "split into", count[ln], "non-overlapping regions."
}' "$OUTDIR/gtf2fasta_nonoverlap_vs_gencode.bed" | tee "$OUTDIR/split_cds_summary.txt"

# === Step 6: Expand all non-overlapping portions into independent GTF entries ===
echo "Step 6: Creating expanded GTF (one entry per non-overlapping portion)..."

awk 'BEGIN{OFS="\t"}
NR==FNR { gtf_line[NR]=$0; next }
{
    original_line_num=$4
    new_chr=$1
    new_start=$2+1
    new_end=$3
    new_strand=$6
    split(gtf_line[original_line_num], fields, "\t")
    print new_chr, "GTF2FASTA", "CDS", new_start, new_end, ".", new_strand, ".", fields[9]
}' "$OUTDIR/gtf2fasta_cds.gtf" "$OUTDIR/gtf2fasta_nonoverlap_vs_gencode.bed" \
> "$OUTDIR/gtf2fasta_unique_segments_expanded.gtf"

# === Summary ===
echo "📊 Final summary:"
echo "  Original GTF2FASTA CDS entries: $(wc -l < "$OUTDIR/gtf2fasta_cds.gtf")"
echo "  GENCODE CDS (with ±35nt flank): $(wc -l < "$OUTDIR/GENCODEv36_cds_flanked.bed")"
echo "  Non-overlapping BED regions:    $(wc -l < "$OUTDIR/gtf2fasta_nonoverlap_vs_gencode.bed")"
echo "  Unique GTF entries:             $(wc -l < "$OUTDIR/gtf2fasta_unique_segments.gtf")"
echo "  Expanded GTF entries:           $(wc -l < "$OUTDIR/gtf2fasta_unique_segments_expanded.gtf")"

echo -e "\n✅ Output files:"
echo "  • GENCODE flanked BED:  $OUTDIR/GENCODEv36_cds_flanked.bed"
echo "  • BED regions:          $OUTDIR/gtf2fasta_nonoverlap_vs_gencode.bed"
echo "  • Unique GTF entries:   $OUTDIR/gtf2fasta_unique_segments.gtf"
echo "  • Expanded GTF:         $OUTDIR/gtf2fasta_unique_segments_expanded.gtf"
echo "  • Split summary:        $OUTDIR/split_cds_summary.txt"

echo -e "\n🔍 Preview (first 5 lines of expanded GTF):"
head -5 "$OUTDIR/gtf2fasta_unique_segments_expanded.gtf"

# === Check flanking for ENSMUST00000208660.2 ===
echo -e "\n🔍 Original GENCODE CDS for ENSMUST00000208660.2:"
grep 'ENSMUST00000208660.2' "$OUTDIR/GENCODEv36_cds.gtf" | awk -F'\t' '{print "  " $1 ":" $4 "-" $5 " (" $7 ")"}'

echo -e "\n🔍 Flanked BED for ENSMUST00000208660.2:"
grep 'ENSMUST00000208660.2' "$OUTDIR/GENCODEv36_cds_flanked.bed" | awk -F'\t' '{print "  " $1 ":" $2+1 "-" $3 " (" $6 ")"}'

#===============================================================================
# EXAMPLE: Overlapping smORF handling with ±35nt flanking
#===============================================================================
# Original smORF CDS:             |-smORF-|      |-smORF-|
# GENCODE CDS overlap:               |---GENCODE---|
# GENCODE with ±35nt:              |----GENCODE+35----|
# Resulting unique portions:     |-|                  |--|
#
# The smORF is split to exclude the GENCODE CDS plus 35nt flanking regions,
# ensuring ribosome footprints near canonical start/stop codons are not
# double-counted as novel smORF reads.
#===============================================================================
