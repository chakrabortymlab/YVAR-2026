#!/bin/bash
#SBATCH --job-name=map_stellate_x
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8         
#SBATCH --mem=32G                  
#SBATCH --time=01:00:00           
#SBATCH --output=map_ste_x_%j.out
#SBATCH --error=map_ste_x_%j.err

set -euo pipefail
# ==========================================
# 1. SETUP AND MODULES
# ==========================================
module purge
module load GCC/13.2.0 OpenMPI/4.1.6 BLAST+/2.16.0 BEDTools/2.31.1 

STRAINS=("ISO1" "A3" "A4")
STE_QUERIES="HetSte.fasta EuSte.fasta"
COMBINED_QUERY="All_Ste_Combined.fasta"
cat $STE_QUERIES > $COMBINED_QUERY

SUMMARY_LOG="Stellate_X_BLAST_Counts_Fixed.txt"
echo -e "Strain\tStellate_Unit_Count" > $SUMMARY_LOG

for STRAIN in "${STRAINS[@]}"; do
    
    X_FASTA="${STRAIN}HarshX.fasta"
    echo "---------------------------------------------------"
    echo "Processing Strain: $STRAIN"

    if [[ ! -f "$X_FASTA" ]]; then continue; fi

    # Step A: Build BLAST Database
    makeblastdb -in $X_FASTA -dbtype nucl -out ${STRAIN}_X_db

    # Step B: Run BLAST
    # outfmt 6 args: 1:target, 2:start, 3:end, 4:query, 5:identity, 6:length
    echo "  Running BLASTn..."
    blastn -query $COMBINED_QUERY \
           -db ${STRAIN}_X_db \
           -outfmt "6 sseqid sstart send qseqid pident length" \
           -dust no \
           -num_threads $SLURM_CPUS_PER_TASK \
           -out ${STRAIN}_Ste_hits_raw.tsv

    # Step C: Format for BED6, flip strands if necessary, and sort by Identity
    echo "  Formatting hits..."
    awk '$5 >= 90 && $6 >= 800 {
            if ($2 < $3) { print $1"\t"$2-1"\t"$3"\t"$4"\t"$5"\t+" }
            else { print $1"\t"$3-1"\t"$2"\t"$4"\t"$5"\t-" }
         }' ${STRAIN}_Ste_hits_raw.tsv | sort -k5,5nr > ${STRAIN}_sorted_by_identity.bed

    # Step D: The Custom Tandem Array Filter (Discards overlapping redundant hits, keeps adjacent units)
    echo "  Filtering tandem arrays..."
    awk -F'\t' '{
        overlap_found = 0;
        for (i=1; i<=count; i++) {
            if ($1 == chrom[i]) {
                start_max = ($2 > start[i]) ? $2 : start[i];
                end_min = ($3 < end[i]) ? $3 : end[i];
                if (end_min > start_max) {
                    ovl_len = end_min - start_max;
                    # If overlap is > 100bp, discard
                    if (ovl_len > 25) {
                        overlap_found = 1;
                        break;
                    }
                }
            }
        }
        if (overlap_found == 0) {
            count++;
            chrom[count] = $1;
            start[count] = $2;
            end[count] = $3;
            print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6;
        }
    }' ${STRAIN}_sorted_by_identity.bed | sort -k1,1 -k2,2n > ${STRAIN}_Ste_best_footprints.bed

    # Step E: Count 
    # Column 4 contains the best-match query ID (EuSte or HetSte)
    UNIT_COUNT=$(wc -l < ${STRAIN}_Ste_best_footprints.bed)
    EUSTE_COUNT=$(awk '$4 == "EuSte" {print}' ${STRAIN}_Ste_best_footprints.bed | wc -l)
    HETSTE_COUNT=$(awk '$4 == "HetSte" {print}' ${STRAIN}_Ste_best_footprints.bed | wc -l)
    
    echo -e "${STRAIN}\t${UNIT_COUNT}\t${EUSTE_COUNT}\t${HETSTE_COUNT}" >> $SUMMARY_LOG
    echo "  Identified $UNIT_COUNT total units ($EUSTE_COUNT EuSte, $HETSTE_COUNT HetSte)."
    # Step F: Extract Sequences (Added -name+ to carry over EuSte/HetSte ID)
    echo "  Extracting FASTA sequences..."
    bedtools getfasta -fi $X_FASTA \
                      -bed ${STRAIN}_Ste_best_footprints.bed \
                      -s \
                      -name+ \
                      -fo ${STRAIN}_Ste_extracted_raw.fasta

    # Step G: Clean Headers to include Strain, Number, Type, and Coordinates
    echo "  Formatting FASTA headers..."
    awk -v prefix="${STRAIN}_X_Ste_" '/^>/ {
        i++
        # Strip the > character
        raw = substr($0, 2)
        # Split by :: (bedtools name+ delimiter)
        split(raw, parts, "::")
        ste_type = parts[1]
        # Split by ( to remove the (+) or (-) strand markup
        split(parts[2], coords, "(")
        
        # Output clean header: e.g., >A4_X_Ste_1_EuSte_scaffold1:1975002-1977114
        print ">" prefix i "_" ste_type "_" coords[1]
        next
    } {print}' ${STRAIN}_Ste_extracted_raw.fasta > Final_${STRAIN}_X_Stellate_BLAST_Units.fasta

    # Cleanup
    rm ${STRAIN}_X_db* ${STRAIN}_Ste_hits_raw.tsv ${STRAIN}_sorted_by_identity.bed ${STRAIN}_Ste_extracted_raw.fasta

done

echo "---------------------------------------------------"
cat $SUMMARY_LOG
date
