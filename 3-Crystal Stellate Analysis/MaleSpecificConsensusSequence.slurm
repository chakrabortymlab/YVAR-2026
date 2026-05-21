#!/bin/bash
#SBATCH --job-name=dRNA_profiler
#SBATCH --output=dRNA_%j.out
#SBATCH --error=dRNA_%j.err
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=168G
#SBATCH --time=06:00:00
#SBATCH --output=dRNA_profiler_%j.log

set -euo pipefail

module purge
module load GCC/13.2.0 minimap2/2.29 SAMtools/1.21 Pysam/0.22.0 SciPy-bundle/2023.11
# ==============================================================================
# 2. DEFINE INPUTS & REFERENCES
# ==============================================================================
MALE_RNA="/scratch/group/chakraborty_lab/adolfo/fertilitygenes/TranscriptData/A4_M_wholebodyRNA_ONT.fq.gz"
FEMALE_RNA="/scratch/group/chakraborty_lab/adolfo/fertilitygenes/TranscriptData/A4_F_wholebodyRNA_ONT.fq.gz"

PCKR_REF="ISO-1PCKR.fasta"
SUSTE_REF="ISO-1SuSte.fasta"

REPORT_FILE="ISO1A4_Transcriptomics_Full_Report.txt"

# ==============================================================================
# 3. GENERATE THE PYSAM PYTHON SCRIPT ON THE FLY
# ==============================================================================
cat << 'EOF' > analyze_bam.py
import sys
import pysam
import numpy as np
from collections import Counter

bam_file = sys.argv[1]
report_file = sys.argv[2]
prefix = sys.argv[3]

# Open BAM File
bam = pysam.AlignmentFile(bam_file, "rb")

total_reads = 0
sense_reads = 0
antisense_reads = 0
read_lengths = []
introns = []

for read in bam.fetch():
    total_reads += 1
    read_lengths.append(read.query_alignment_length)
    
    # Check strand orientation
    if read.is_reverse:
        antisense_reads += 1
    else:
        sense_reads += 1
        
    # Parse CIGAR for introns (Operation 3 is 'N')
    current_ref_pos = read.reference_start
    for operation, length in read.cigartuples:
        if operation == 3: 
            if length > 30: # Filter out small sequencing indels
                introns.append((current_ref_pos, current_ref_pos + length, length))
            current_ref_pos += length
        elif operation in (0, 2, 7, 8): # M, D, =, X (consume reference)
            current_ref_pos += length

bam.close()

# Write output to the shared report file
with open(report_file, "a") as f:
    f.write(f"===========================================================\n")
    f.write(f"Processing: {prefix}\n")
    f.write(f"===========================================================\n")
    f.write(f"Total Primary Mapped Reads: {total_reads}\n")
    
    if total_reads == 0:
        f.write("No strictly mapped reads found in this library.\n\n")
        sys.exit(0)
        
    f.write(f"Strand Bias:\n")
    f.write(f"  - Sense (+):     {sense_reads}\n")
    f.write(f"  - Antisense (-): {antisense_reads}\n")
    
    f.write(f"Read Lengths:\n")
    f.write(f"  - Max Length:  {max(read_lengths)} bp\n")
    f.write(f"  - Mean Length: {int(np.mean(read_lengths))} bp\n")
    
    f.write(f"Structural Splicing (Intron Architecture):\n")
    if not introns:
        f.write("  - No introns detected. Appears to be a single-exon transcript (or unspliced).\n")
    else:
        # Group identical introns and count supporting reads
        intron_counts = Counter(introns)
        for (start, end, length), count in sorted(intron_counts.items()):
            f.write(f"  - Span {start} to {end} (Length: {length} bp) [Supported by {count} reads]\n")
    f.write("\n")
EOF

# ==============================================================================
# 4. MAIN MAPPING & PROFILING FUNCTION
# ==============================================================================
map_and_profile() {
    local ref=$1
    local reads=$2
    local prefix=$3

    echo "Aligning $prefix..."
    
    # 1. Map with minimap2
    minimap2 -t 16 -ax splice -uf -k14 "$ref" "$reads" > "${prefix}_raw.sam" 

    # 2. Strict Filter & Sort with samtools (Exclude unmapped, secondary, supplementary)
    samtools view -@ 16 -bS -F 2308 -q 60 "${prefix}_raw.sam" | samtools sort -@ 16 -o "${prefix}_strict.bam"
    samtools index -@ 16 "${prefix}_strict.bam"
    
    # 3. Generate base-by-base coverage map for later plotting
    samtools depth -a "${prefix}_strict.bam" > "${prefix}_coverage.txt"
    rm "${prefix}_raw.sam"

    # 4. Run the Python Profiler
    python analyze_bam.py "${prefix}_strict.bam" "$REPORT_FILE" "$prefix"
}

# ==============================================================================
# 5. EXECUTE PIPELINE
# ==============================================================================
echo "Starting dRNA-seq Pipeline with Pysam..."

map_and_profile "$PCKR_REF" "$MALE_RNA" "ISO1A4_Male_PCKR"
map_and_profile "$PCKR_REF" "$FEMALE_RNA" "ISO1A4_Female_PCKR"
map_and_profile "$SUSTE_REF" "$MALE_RNA" "ISO1A4_Male_SuSte"
map_and_profile "$SUSTE_REF" "$FEMALE_RNA" "ISO1A4_Female_SuSte"

echo "Pipeline complete. All details saved to $REPORT_FILE."
