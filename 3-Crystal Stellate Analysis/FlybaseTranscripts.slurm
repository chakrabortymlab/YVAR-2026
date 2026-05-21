#!/usr/bin/env bash
#SBATCH --job-name=flybase_blast
#SBATCH --output=flybase_blast.log
#SBATCH --time=05:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=248G

set -euo pipefail

# ---1. Load BLAST Module ---
module purge
module load GCC/13.3.0  OpenMPI/5.0.3 BLAST+/2.16.0

echo "=== Starting FlyBase Validation Pipeline ==="

FB_FASTA="dmel-all-transcript-r6.67.fasta"

# --- 2. Build BLAST Database ---
echo "Building local BLAST database..."
makeblastdb -in "$FB_FASTA" -dbtype nucl -out fb_transcript_db

# --- 3. Prepare Queries ---
# Combine your consensus sequences into a single query file.
# Assumes you have PCKR_consensus.fasta, Ssl.fasta, and SuSte.fasta in the directory.
cat ISO-1PCKR.fasta ISO-1SuSte.fasta HetSte.fasta EuSte.fasta Ssl.fasta > all_queries.fasta

# --- 4. Run High-Stringency BLASTn ---
echo "Running BLASTn against FlyBase..."
# -outfmt 6 gives a clean tabular output. 
# We add the sequence titles (stitle) so we can see the actual FlyBase gene names (e.g., CG40635).
blastn -query all_queries.fasta \
       -db fb_transcript_db \
       -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
       -evalue 1e-10 \
       -num_threads 32 \
       -max_target_seqs 5 \
       -out FlyBase_Validation_ResultsTranscriptAllISO1.tsv

echo "=== BLAST Complete ==="
echo "Top hits for PCKR:"
grep -i "PCKR" FlyBase_Validation_Results.tsv | head -n 5

echo "Top hits for Ssl:"
grep -i "Ssl" FlyBase_Validation_Results.tsv | head -n 5

echo "Top hits for SuSte:"
grep -i "SuSte" FlyBase_Validation_Results.tsv | head -n 5
