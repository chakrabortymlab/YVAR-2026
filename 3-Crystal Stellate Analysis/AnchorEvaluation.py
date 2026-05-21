import os
import glob
import re
import pandas as pd
from itertools import combinations

# ==========================================
# CONFIGURATION
# ==========================================
ALN_DIR = "aligned_anchors"
OUTPUT_TSV = "Conserved_Anchors_Locations.tsv"
OUTPUT_SUMMARY = "Conserved_Inverted_Anchors_Summary.csv"

# Minimum Average Pairwise Identity to be considered a "Conserved Anchor"
MIN_IDENTITY = 95.0  

def parse_fasta(filepath):
    """Parses a FASTA alignment file."""
    seqs = {}
    with open(filepath, 'r') as f:
        header = None
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                header = line
                seqs[header] = []
            else:
                seqs[header].append(line)
    return {k: "".join(v) for k, v in seqs.items()}

def calc_pairwise_identity(seq1, seq2):
    """Calculates identity between two aligned sequences, ignoring double-gaps."""
    matches = 0
    valid_length = 0
    for a, b in zip(seq1, seq2):
        if a == '-' and b == '-': 
            continue
        valid_length += 1
        if a.upper() == b.upper(): 
            matches += 1
    return (matches / valid_length) * 100 if valid_length > 0 else 0.0

def evaluate_alignments():
    aln_files = glob.glob(os.path.join(ALN_DIR, "*.aln"))
    print(f"Found {len(aln_files)} alignment files to evaluate...")
    
    valid_anchors = []
    summary_data = []
    
    for aln in aln_files:
        group_id = os.path.basename(aln).replace(".aln", "")
        seqs = parse_fasta(aln)
        
        if len(seqs) != 3:
            continue # Skip if it doesn't have exactly 3 strains
            
        # Calculate pairwise identities
        identities = []
        for (h1, s1), (h2, s2) in combinations(seqs.items(), 2):
            identities.append(calc_pairwise_identity(s1, s2))
            
        avg_identity = sum(identities) / len(identities)
        
        # If the trio is highly conserved, extract coordinates and save
        if avg_identity >= MIN_IDENTITY:
            summary_data.append({'Group': group_id, 'Avg_Identity': round(avg_identity, 2)})
            
            for header in seqs.keys():
                # Extract coordinates from header, ignoring optional strand flag like (-) or (+)
                match = re.search(r'>([A-Z0-9\-]+)_([^:]+):(\d+)-(\d+)(?:\([+-]\))?_RelPos_([0-9\.]+)', header)
                if match:
                    valid_anchors.append({
                        'Group_ID': group_id,
                        'Strain': match.group(1),
                        'Scaffold': match.group(2),
                        'Anchor_Start': int(match.group(3)),
                        'Anchor_End': int(match.group(4)),
                        'Avg_Identity': round(avg_identity, 2),
                        'Relative_Pos': float(match.group(5))
                    })
                else:
                    print(f"Failed to parse header: {header}")

    if valid_anchors:
        df_tsv = pd.DataFrame(valid_anchors)
        df_tsv = df_tsv.sort_values(by=['Group_ID', 'Strain'])
        df_tsv.to_csv(OUTPUT_TSV, sep='\t', index=False)
        
        df_summary = pd.DataFrame(summary_data)
        df_summary.to_csv(OUTPUT_SUMMARY, index=False)
        
        print(f"\nSuccess! Kept {len(df_summary)} highly conserved inverted anchor groups (>={MIN_IDENTITY}% identity).")
        print(f"Coordinate TSV saved to: {OUTPUT_TSV}")
        print(f"Summary CSV saved to: {OUTPUT_SUMMARY}")
    else:
        print(f"\nNo groups passed the strict {MIN_IDENTITY}% sequence identity threshold.")

if __name__ == "__main__":
    evaluate_alignments()
