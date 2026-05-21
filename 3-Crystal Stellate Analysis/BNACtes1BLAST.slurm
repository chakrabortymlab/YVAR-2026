#!/bin/bash
#SBATCH --job-name=flybase_genes
#SBATCH --output=flybase_genes_%j.out
#SBATCH --time=00:30:00             
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16           
#SBATCH --mem=64G                    

module purge
module load GCC/13.3.0  OpenMPI/5.0.3 BLAST+/2.16.0

echo "Building FlyBase Extended Genes Database..."
makeblastdb -in dmel-all-gene_extended2000-r6.67.fasta -dbtype nucl -out flybase_extended_db

OUTFMT="6 qseqid sseqid pident length qstart qend sstart send evalue stitle"

echo "1. Testing ISO-1 Su(Ste)..."
blastn -query ISO-1SuSte.fasta \
       -db flybase_extended_db \
       -out SuSte_vs_Genes.tsv \
       -outfmt "$OUTFMT" \
       -evalue 1e-5 \
       -word_size 11 \
       -num_threads $SLURM_CPUS_PER_TASK

echo "2. Testing ISO-1 PCKR..."
blastn -query ISO-1PCKR.fasta \
       -db flybase_extended_db \
       -out PCKR_vs_Genes.tsv \
       -outfmt "$OUTFMT" \
       -evalue 1e-5 \
       -word_size 11 \
       -num_threads $SLURM_CPUS_PER_TASK

echo "3. Testing Heterochromatic Stellate (HetSte)..."
blastn -query HetSte.fasta \
       -db flybase_extended_db \
       -out HetSte_vs_Genes.tsv \
       -outfmt "$OUTFMT" \
       -evalue 1e-5 \
       -word_size 11 \
       -num_threads $SLURM_CPUS_PER_TASK

echo "4. Testing Euchromatic Stellate (EuSte)..."
blastn -query EuSte.fasta \
       -db flybase_extended_db \
       -out EuSte_vs_Genes.tsv \
       -outfmt "$OUTFMT" \
       -evalue 1e-5 \
       -word_size 11 \
       -num_threads $SLURM_CPUS_PER_TASK

echo "All BLAST searches complete!"

echo "Filtering for target genes (βNACtes1, Ssl, etc.)..."

# Create a combined master file for ALL filtered hits across the 4 searches
MASTER_FILTERED="All_Filtered_Targets.tsv"
> "$MASTER_FILTERED" # Initialize empty file

# Loop through the 4 output TSVs to filter them
for FILE in SuSte_vs_Genes.tsv PCKR_vs_Genes.tsv HetSte_vs_Genes.tsv EuSte_vs_Genes.tsv; do
    
    # Define the output filename for the specific filtered query
    PREFIX="${FILE%_vs_Genes.tsv}"
    FILTERED_FILE="${PREFIX}_Filtered.tsv"
    
    # -i makes it case-insensitive
    # -E allows extended regular expressions (the | OR operator)
    grep -iE "FBgn0030538|βNACtes1|CG13402|FBgn0015300|Ssl|CG13591" "$FILE" > "$FILTERED_FILE"
    
    # Append these hits to the master file as well
    cat "$FILTERED_FILE" >> "$MASTER_FILTERED"
    
    # Count how many hits were found and print it to the SLURM output log
    HIT_COUNT=$(wc -l < "$FILTERED_FILE")
    echo " -> Found $HIT_COUNT target hits in $FILE"
done

echo "Done! Target hits saved to individual _Filtered.tsv files and a combined $MASTER_FILTERED file."
