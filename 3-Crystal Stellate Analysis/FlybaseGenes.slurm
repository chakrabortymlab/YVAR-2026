#!/bin/bash
#SBATCH --job-name=flybase_genes
#SBATCH --output=flybase_genes_%j.out
#SBATCH --time=01:30:00             
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16           
#SBATCH --mem=64G

set -euo pipefail

module purge
module load GCC/13.3.0  OpenMPI/5.0.3 BLAST+/2.16.0

# The query file containing all genes + 2kb flanking sequence to capture promoters
QUERY_FASTA="dmel-all-gene_extended2000-r6.67.fasta"

# BLAST output format
OUTFMT="6 qseqid sseqid pident length qstart qend sstart send evalue stitle"

for PREFIX in "A4" "A3" "ISO-1"; do
    ARRAY_FASTA="${PREFIX}_Array.fasta"
    ORIG_BED="${PREFIX}.all_hits.bed"

    if [[ ! -f "$ARRAY_FASTA" || ! -f "$ORIG_BED" ]]; then
        echo "WARNING: Missing files for $PREFIX. Skipping..."
        continue
    fi

    echo "========================================"
    echo "Processing $PREFIX Array..."
    echo "========================================"

    # 1. Extract the exact FASTA header to use as the IGV Chromosome Name
    IGV_CHROM=$(head -n 1 "$ARRAY_FASTA" | sed 's/>//')
    
    # Extract the starting coordinate to use as the offset
    OFFSET=$(echo "$IGV_CHROM" | awk -F'[:-]' '{print $2}')

    echo "Detected IGV Chromosome: $IGV_CHROM"
    echo "Coordinate Offset: $OFFSET"

    # 2. Shift the original BED annotations down to the Array's coordinate space
    awk -v chr="$IGV_CHROM" -v ofs="$OFFSET" -F'\t' 'BEGIN{OFS="\t"} {
        if ($1 != "" && $2 ~ /^[0-9]+$/) {
            new_start = $2 - ofs
            new_end = $3 - ofs
            if (new_start < 0) new_start = 0
            print chr, new_start, new_end, $4, $5, $6
        }
    }' "$ORIG_BED" > "${PREFIX}_Array_Annotations.bed"

    echo "Generated shifted annotation BED: ${PREFIX}_Array_Annotations.bed"

    # 3. Build BLAST Database from the Array (FIXED: removed -quiet, redirected output)
    DB_NAME="${PREFIX}_db"
    makeblastdb -in "$ARRAY_FASTA" -dbtype nucl -out "$DB_NAME" > /dev/null 2>&1

    # 4. Run BLAST search
    TSV_OUT="${PREFIX}_All_Genes.tsv"
    echo "Running BLAST for $PREFIX..."
    blastn -query "$QUERY_FASTA" -db "$DB_NAME" -out "$TSV_OUT" -outfmt "$OUTFMT" -evalue 1e-10 -num_threads $SLURM_CPUS_PER_TASK

    # 5. Extract ALL hits into a master BED file
    awk -v chr="$IGV_CHROM" -F'\t' 'BEGIN{OFS="\t"} {
        sstart=$7; send=$8; pident=$3; qseqid=$1
        if (sstart < send) {
            start = sstart - 1
            end = send
            strand = "+"
        } else {
            start = send - 1
            end = sstart
            strand = "-"
        }
        name = qseqid "_id:" pident "%"
        print chr, start, end, name, ".", strand
    }' "$TSV_OUT" > "${PREFIX}_Array_All_BLAST.bed"
    
    echo "Generated master BLAST BED: ${PREFIX}_Array_All_BLAST.bed"

    # 6. Extract the specific genes of interest (BNACtes and Ssl)
    grep -E "FBgn0030538|βNACtes1|CG13402" "$TSV_OUT" > "${PREFIX}_BNACtes.tsv"
    grep -E "FBgn0015300|Ssl|CG13591" "$TSV_OUT" > "${PREFIX}_Ssl.tsv"

    # 7. Convert filtered TSVs to BED
    awk -v chr="$IGV_CHROM" -F'\t' 'BEGIN{OFS="\t"} {
        sstart=$7; send=$8; pident=$3; qseqid=$1
        if (sstart < send) {
            start = sstart - 1
            end = send
            strand = "+"
        } else {
            start = send - 1
            end = sstart
            strand = "-"
        }
        name = qseqid "_id:" pident "%"
        print chr, start, end, name, ".", strand
    }' "${PREFIX}_BNACtes.tsv" > "${PREFIX}_Array_BNACtes_Promoters.bed"

    awk -v chr="$IGV_CHROM" -F'\t' 'BEGIN{OFS="\t"} {
        sstart=$7; send=$8; pident=$3; qseqid=$1
        if (sstart < send) {
            start = sstart - 1
            end = send
            strand = "+"
        } else {
            start = send - 1
            end = sstart
            strand = "-"
        }
        name = "Ssl_id:" pident "%"
        print chr, start, end, name, ".", strand
    }' "${PREFIX}_Ssl.tsv" > "${PREFIX}_Array_Ssl_Promoters.bed"

    # Cleanup intermediate files
    rm "${PREFIX}_BNACtes.tsv" "${PREFIX}_Ssl.tsv"

    echo "Finished $PREFIX. BED files ready for IGV."
done

echo "All processing complete!"
