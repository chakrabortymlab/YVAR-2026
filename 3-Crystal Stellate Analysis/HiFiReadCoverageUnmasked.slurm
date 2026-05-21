#!/bin/bash
#SBATCH --job-name=winnowmap_array        
#SBATCH --ntasks-per-node=1                          	
#SBATCH --cpus-per-task=48                 
#SBATCH --mem=164G                          
#SBATCH --time=24:00:00                   
#SBATCH --output=winnowmap_pipeline_%j.out 
#SBATCH --error=winnowmap_pipeline_%j.err  

set -euo pipefail

# --- Configuration ---
declare -A READS
READS=(
    ["ISO1"]="/scratch/user/aadelgado6/Delgado/fertilitygenes/readconfirmation4codingsequences/ISO1HiFireads.fastq.gz"
    ["A3"]="/scratch/user/aadelgado6/Delgado/fertilitygenes/readconfirmation4codingsequences/A3HiFireads.fastq.gz"
    ["A4"]="/scratch/user/aadelgado6/Delgado/fertilitygenes/readconfirmation4codingsequences//A4HiFireads.fastq.gz"
)

# Dynamically inherit threads from Slurm allocation
THREADS=48
KMER_SIZE=21
HASH_SIZE="8G" # Memory allocation for Jellyfish

echo "Starting Winnowmap2 Array Pipeline on $(hostname) at $(date)"

for STRAIN in "${!READS[@]}"; do
    echo "========================================================="
    echo " Processing strain: $STRAIN"
    echo "========================================================="

    REF="${STRAIN}_Array.fasta"
    HIFI_READS="${READS[$STRAIN]}"
    GAP_REPORT="${STRAIN}_gap_report.tsv"
    MERYL_DB="${STRAIN}_merylDB"
    REP_KMERS="${STRAIN}_repetitive_k15.txt"
    BAM_OUT="${STRAIN}_array_strict.bam"
    FQ_OUT="${STRAIN}_array_reads.fastq"
    JF_OUT="${STRAIN}_k${KMER_SIZE}.jf"
    HISTO_OUT="${STRAIN}_k${KMER_SIZE}.histo"

    # Step 1: Identify Gaps (Ns) in the Array using pure awk (HPC optimized)
    echo "[1/5] Scanning $REF for assembly gaps (Ns)..."
    echo -e "Contig\tGap_Start\tGap_End\tGap_Size(bp)" > "$GAP_REPORT"
    
    awk '
    /^>/ {
        if (name != "") {
            offset = 0;
            temp_seq = seq;
            while(match(temp_seq, /N+/)) {
                start = offset + RSTART;
                len = RLENGTH;
                end = start + len - 1;
                printf "%s\t%d\t%d\t%d\n", name, start, end, len;
                temp_seq = substr(temp_seq, RSTART + RLENGTH);
                offset = end;
            }
        }
        name = substr($1, 2); 
        seq = "";             
        next;
    }
    {
        seq = seq toupper($0);
    }
    END {
        if (name != "") {
            offset = 0;
            temp_seq = seq;
            while(match(temp_seq, /N+/)) {
                start = offset + RSTART;
                len = RLENGTH;
                end = start + len - 1;
                printf "%s\t%d\t%d\t%d\n", name, start, end, len;
                temp_seq = substr(temp_seq, RSTART + RLENGTH);
                offset = end;
            }
        }
    }' "$REF" >> "$GAP_REPORT"
    
    GAP_COUNT=$(tail -n +2 "$GAP_REPORT" | wc -l)
    echo "      -> Found $GAP_COUNT gaps. Details saved to $GAP_REPORT."

    # Step 2: Build the Meryl k-mer database for Winnowmap
    module purge
    module load GCC/13.2.0 Meryl/1.4.1
    echo "[2/5] Building Meryl database for repetitive k-mers..." 
    meryl count k=15 output "$MERYL_DB" "$REF"
    meryl print greater-than distinct=0.9998 "$MERYL_DB" > "$REP_KMERS"

    # Step 3: Map reads strictly using Winnowmap2
    echo "[3/5] Mapping HiFi reads with Winnowmap2..."
    module purge
    module load GCC/12.3.0 Winnowmap/2.03 SAMtools/1.18
    winnowmap -W "$REP_KMERS" -t "$THREADS" -ax map-pb --secondary=no "$REF" "$HIFI_READS" | \
        samtools view -@ "$THREADS" -bS -F 2308 - | \
        samtools sort -@ "$THREADS" -o "$BAM_OUT"
    
    samtools index "$BAM_OUT"

    # Step 4: Extract the mapped reads back to FASTQ format
    echo "[4/5] Extracting mapped reads..."
    samtools fastq -@ "$THREADS" "$BAM_OUT" > "$FQ_OUT"

    # Step 5: Run Jellyfish to count canonical k-mers
    echo "[5/5] Counting k-mers with Jellyfish & generating GenomeScope histo..."
    module purge
    module load GCC/13.2.0 Jellyfish/2.3.1
    jellyfish count -C -m "$KMER_SIZE" -s "$HASH_SIZE" -t "$THREADS" -o "$JF_OUT" "$FQ_OUT"
    jellyfish histo -t "$THREADS" "$JF_OUT" > "$HISTO_OUT"

    echo ">>> $STRAIN complete."
    echo "---------------------------------------------------------"
done

echo "Pipeline complete at $(date)!"
