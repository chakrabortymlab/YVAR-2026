#!/usr/bin/env bash
#SBATCH --job-name=dual_array
#SBATCH --output=dual_array.%j.log
#SBATCH --time=1:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G

set -euo pipefail

########################################
### 1. USER INPUTS
########################################
STRAIN_ID="A3"
ASSEMBLY="DelgadoYA3.fasta"
READS_BAM="mappedreads.bam"
RM_OUT="DelgadoYA3.fasta.out" # RepeatMasker output for flanking analysis

SUSTE_SEQ="SuSte.fasta"
PCKR_SEQ="PCKR_consensus.fasta"

# Empirically Determined Family-Specific Coverages from IGV read alignment visualization
COV_SUSTE=42
COV_PCKR=33

# Expected unit lengths for CNV math
LEN_SUSTE=2806
LEN_PCKR=2563

module purge
module load GCC/13.2.0 SAMtools/1.21 minimap2/2.29 BEDTools/2.31.1

echo "=== Starting Su(Ste) & PCKR Dual Analysis for ${STRAIN_ID} ==="

cat "$SUSTE_SEQ" "$PCKR_SEQ" > dual_queries.fasta
samtools faidx "$ASSEMBLY"

########################################
### 2. MAPPING & FILTERING
########################################
echo "[1/4] Mapping Consensuses (Repetitive Filters Disabled)..."
# -f 0 and -P force minimap2 to evaluate ALL repetitive k-mers and retain all chains
minimap2 -t "${SLURM_CPUS_PER_TASK:-8}" -f 0 -P -N 10000 -p 0.5 -c "${ASSEMBLY}" dual_queries.fasta \
    > "${STRAIN_ID}.raw.paf"

echo "[2/4] Filtering Hits (Target Length >= 1850, Query-based PID >= 80)..."
awk 'BEGIN{OFS="\t"} {
    q_len = $4 - $3;
    t_len = $9 - $8;
    
    # Calculate identity based ONLY on the mapped portion of the consensus query.
    # This perfectly tolerates large TE insertions in the target genome.
    pid = ($10 / q_len) * 100;  
    
    if (t_len >= 1850 && pid >= 80) {
        print $6, $8, $9, $1, t_len, $5, pid
    }
}' "${STRAIN_ID}.raw.paf" \
    | sort -k1,1 -k2,2n \
    > "${STRAIN_ID}.all_hits.bed"

# Split arrays
grep 'SuSte' "${STRAIN_ID}.all_hits.bed" > "${STRAIN_ID}.suste.bed" || true
grep 'PCKR'  "${STRAIN_ID}.all_hits.bed" > "${STRAIN_ID}.pckr.bed" || true

########################################
### 3. COPY NUMBER & ARCHITECTURE
########################################
echo "[3/4] Calculating Copy Number..."
REPORT_FILE="${STRAIN_ID}_DUAL_REPORT.txt"
cat <<EOF > "$REPORT_FILE"
########################################################
### SU(STE) / PCKR DUAL ARRAY REPORT: ${STRAIN_ID}
########################################################
EOF

analyze_family() {
    local family=$1; local bed=$2; local len=$3; local depth=$4
    if [[ ! -s "$bed" ]]; then return; fi
    
    local vis=$(wc -l < "$bed")
    local avg_id=$(awk '{sum += $7} END {printf "%.2f", sum/NR}' "$bed")
    local plus=$(awk '$6 == "+"' "$bed" | wc -l)
    local minus=$(awk '$6 == "-"' "$bed" | wc -l)
    
    bedtools merge -i "$bed" > "${bed}.merged"
    samtools bedcov "${bed}.merged" "${READS_BAM}" > "${bed}.mass"
    local total_mass=$(awk '{sum += $NF} END {print sum}' "${bed}.mass")
    
    # Calculate estimates using the family-specific coverage variable ($depth)
    local est=$(awk -v m="$total_mass" -v d="$depth" -v l="$len" \
        'BEGIN { if (d > 0) printf "%.2f", m/(d*l); else print "0" }')
    local collapse=$(awk -v e="$est" -v v="$vis" \
        'BEGIN { if (v > 0) printf "%.2f", e/v; else print "0" }')

    cat <<EOF >> "$REPORT_FILE"
--------------------------------------------------------
 FAMILY: $family
--------------------------------------------------------
Coverage Baseline: ${depth}x
Visible Units:  $vis
Estimated:      $est
Collapse:       ${collapse}x
Identity:       ${avg_id}%
Strand (+/-):   ${plus}:${minus}
EOF
}

# Pass the specific coverage limits into the analysis function
analyze_family "Su(Ste)" "${STRAIN_ID}.suste.bed" "$LEN_SUSTE" "$COV_SUSTE"
analyze_family "PCKR" "${STRAIN_ID}.pckr.bed" "$LEN_PCKR" "$COV_PCKR"

########################################
### 4. FLANKING REPEAT CONTEXT
########################################
echo "[4/4] Analyzing Flanking Repeats..."

# Merge both arrays into one array footprint
cat "${STRAIN_ID}.suste.bed" "${STRAIN_ID}.pckr.bed" | sort -k1,1 -k2,2n > "${STRAIN_ID}.combined.bed"
bedtools cluster -d 50000 -i "${STRAIN_ID}.combined.bed" > "${STRAIN_ID}.clustered.bed"
bedtools merge -i "${STRAIN_ID}.clustered.bed" > "${STRAIN_ID}.megacluster_footprint.bed"

# Format RM output
awk 'NR>3 && NF>10 {print $5"\t"$6"\t"$7"\t"$11"\t0\t+"}' "${RM_OUT}" \
    | sort -k1,1 -k2,2n > temp_rm.bed

# Get 10kb flanks, subtract the array itself, and intersect with TEs
bedtools slop -i "${STRAIN_ID}.megacluster_footprint.bed" -g "${ASSEMBLY}.fai" -b 10000 \
    | bedtools subtract -a - -b "${STRAIN_ID}.megacluster_footprint.bed" \
    | bedtools intersect -a - -b temp_rm.bed -wa -wb \
    > "${STRAIN_ID}.repeat_context.tsv" || true

awk '{print $8}' "${STRAIN_ID}.repeat_context.tsv" | sort | uniq -c | sort -rn > "${STRAIN_ID}.flanking_repeats.txt" || true
rm -f temp_rm.bed "${STRAIN_ID}.combined.bed" "${STRAIN_ID}.clustered.bed" dual_queries.fasta

cat <<EOF >> "$REPORT_FILE"
--------------------------------------------------------
 FLANKING TEs (Top 5 within 10kb of Array Edges)
--------------------------------------------------------
$(head -n 5 "${STRAIN_ID}.flanking_repeats.txt" || echo "None found")
########################################################
EOF

cat "$REPORT_FILE"
echo "Analysis Complete."
