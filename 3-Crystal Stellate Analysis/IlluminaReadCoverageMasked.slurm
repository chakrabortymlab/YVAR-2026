#!/bin/bash
#SBATCH --job-name=bwa_array_masked_illumina
#SBATCH --ntasks-per-node=1                                          
#SBATCH --cpus-per-task=32                 
#SBATCH --mem=164G                          
#SBATCH --time=10:00:00                   
#SBATCH --output=bwa_array_masked_illumina_%j.out 
#SBATCH --error=bwa_array_masked_illumina_%j.err  

set -euo pipefail

# --- Configuration ---
declare -A READS_R1
READS_R1=(
    ["ISO1"]="./IlluminaReads/ISO1/SRR29479671_1.fastq"
    ["A3"]="./IlluminaReads/A3/SRR29479669_1.fastq"
    ["A4"]="./IlluminaReads/A4/SRR29479670_1.fastq"
)

declare -A READS_R2
READS_R2=(
    ["ISO1"]="./IlluminaReads/ISO1/SRR29479671_2.fastq"
    ["A3"]="./IlluminaReads/A3/SRR29479669_2.fastq"
    ["A4"]="./IlluminaReads/A4/SRR29479670_2.fastq"
)

# Explicitly map to the Hard-Masked versions of the ISOLATED ARRAYS
declare -A REFS
REFS=(
    ["ISO1"]="./RepeatMaskerRun/ISO1_Array.fasta.masked"
    ["A3"]="./RepeatMaskerRun/A3_Array.fasta.masked"
    ["A4"]="./RepeatMaskerRun/A4_Array.fasta.masked"
)

THREADS=32
KMER_SIZE=21
HASH_SIZE="8G" # Memory allocation for Jellyfish

# Create the output directory
OUT_DIR="IlluminaCoverage"
mkdir -p "$OUT_DIR"

echo "Starting Illumina Masked Array Pipeline on $(hostname) at $(date)"
echo "All outputs will be routed to: ./${OUT_DIR}/"

for STRAIN in "${!READS_R1[@]}"; do
    echo "========================================================="
    echo " Processing strain: $STRAIN (Illumina)"
    echo "========================================================="

    REF="${REFS[$STRAIN]}" 
    R1="${READS_R1[$STRAIN]}"
    R2="${READS_R2[$STRAIN]}"
    
    # Output names routed to the IlluminaCoverage directory
    BAM_OUT="${OUT_DIR}/${STRAIN}_illumina_array_masked_strict.bam"
    FQ_OUT="${OUT_DIR}/${STRAIN}_illumina_array_masked_reads.fastq"
    JF_OUT="${OUT_DIR}/${STRAIN}_illumina_array_masked_k${KMER_SIZE}.jf"
    HISTO_OUT="${OUT_DIR}/${STRAIN}_illumina_array_masked_k${KMER_SIZE}.histo"

    if [[ ! -f "$REF" ]]; then
        echo "Error: Cannot find masked array file $REF for $STRAIN! Skipping..."
        continue
    fi

    # Step 1: Index the masked reference for BWA (if not already done)
    if [[ ! -f "${REF}.bwt" ]]; then
        echo "[1/4] Building BWA index for $STRAIN..."
        module purge
        module load bwa-mem2/2.2.1-Linux64
        module load GCC/13.2.0 SAMtools/1.21
        bwa-mem2 index "$REF"
    else
        echo "[1/4] BWA index already exists for $STRAIN."
    fi

    # Step 2: Map reads strictly using BWA-MEM
    echo "[2/4] Mapping Illumina reads with BWA-MEM to Masked Reference..."
    
    # Map, filter for primary alignments (-F 2308), and sort
    bwa-mem2 mem -t "$THREADS" "$REF" "$R1" "$R2" | \
        samtools view -@ "$THREADS" -bS -F 2308 - | \
        samtools sort -@ "$THREADS" -o "$BAM_OUT"
    
    samtools index "$BAM_OUT"

    # Step 3: Extract the mapped reads back to FASTQ format
    echo "[3/4] Extracting mapped reads..."
    samtools fastq -@ "$THREADS" "$BAM_OUT" > "$FQ_OUT"

    # Step 4: Run Jellyfish to count canonical k-mers
    echo "[4/4] Counting k-mers with Jellyfish & generating GenomeScope histo..."
    module purge
    module load GCC/13.2.0 Jellyfish/2.3.1
    jellyfish count -C -m "$KMER_SIZE" -s "$HASH_SIZE" -t "$THREADS" -o "$JF_OUT" "$FQ_OUT"
    jellyfish histo -t "$THREADS" "$JF_OUT" > "$HISTO_OUT"

    echo ">>> $STRAIN complete."
    echo "---------------------------------------------------------"
done

