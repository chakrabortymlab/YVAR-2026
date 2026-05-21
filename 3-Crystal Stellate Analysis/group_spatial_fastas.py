import re
import os
import pandas as pd

# 1. Dynamically calculate TRUE array boundaries from the master BED files
def get_true_bounds(filepath):
    try:
        df = pd.read_csv(filepath, sep='\t', header=None)
        # Handle different bed column lengths safely
        if len(df.columns) < 6:
            df.columns = ['Scaffold', 'Start', 'End', 'Family', 'Score'][:len(df.columns)]
        else:
            df.columns = ['Scaffold', 'Start', 'End', 'Family', 'Score', 'Strand'] + list(df.columns[6:])
            
        search_term = r'PCKR|SuSte|Su\(Ste\)'
        target_df = df[df['Family'].astype(str).str.contains(search_term, case=False, na=False, regex=True)].copy()
        
        if target_df.empty: return 0, 0
        
        main_scaffold = target_df['Scaffold'].value_counts().idxmax()
        target_df = target_df[target_df['Scaffold'] == main_scaffold]
        
        return target_df['Start'].min(), target_df['End'].max()
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return 0, 0

iso_min, iso_max = get_true_bounds("ISO-1.all_hits.bed")
a3_min, a3_max = get_true_bounds("A3.all_hits.bed")
a4_min, a4_max = get_true_bounds("A4.all_hits.bed")

def parse_fasta_with_coords(filepath):
    seqs = []
    with open(filepath, 'r') as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    seqs.append({'header': header, 'seq': "".join(seq_lines)})
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            seqs.append({'header': header, 'seq': "".join(seq_lines)})
            
    for s in seqs:
        match = re.search(r':(\d+)-(\d+)', s['header'])
        if match:
            s['start'] = int(match.group(1))
            s['end'] = int(match.group(2))
            s['center'] = (s['start'] + s['end']) / 2
    return seqs

def calculate_relative_positions(seqs, true_min, true_max):
    """Calculates relative position based on the TRUE array boundaries."""
    if not seqs: return []
    span = true_max - true_min
    
    for s in seqs:
        # Calculate depth relative to the strict biological array bounds
        s['rel_pos'] = (s['center'] - true_min) / span if span > 0 else 0.5
    return seqs

# Load and calculate using true boundaries
iso_seqs = calculate_relative_positions(parse_fasta_with_coords("ISO-1_inv.fasta"), iso_min, iso_max)
a3_seqs = calculate_relative_positions(parse_fasta_with_coords("A3_inv.fasta"), a3_min, a3_max)
a4_seqs = calculate_relative_positions(parse_fasta_with_coords("A4_inv.fasta"), a4_min, a4_max)

print(f"Loaded {len(iso_seqs)} ISO-1, {len(a3_seqs)} A3, and {len(a4_seqs)} A4 inverted units.")

# Grouping Logic: Relaxed to 10% to account for lineage-specific array expansions
MAX_SHIFT = 0.10 

group_count = 0
for iso in iso_seqs:
    a3_match = min(a3_seqs, key=lambda x: abs(x['rel_pos'] - iso['rel_pos']), default=None)
    a4_match = min(a4_seqs, key=lambda x: abs(x['rel_pos'] - iso['rel_pos']), default=None)
    
    if a3_match and a4_match:
        dist_a3 = abs(a3_match['rel_pos'] - iso['rel_pos'])
        dist_a4 = abs(a4_match['rel_pos'] - iso['rel_pos'])
        
        if dist_a3 <= MAX_SHIFT and dist_a4 <= MAX_SHIFT:
            group_count += 1
            out_file = f"anchor_groups/inv_anchor_group_{group_count:03d}.fasta"
            with open(out_file, 'w') as out:
                out.write(f">ISO-1_{iso['header']}_RelPos_{iso['rel_pos']:.3f}\n{iso['seq']}\n")
                out.write(f">A3_{a3_match['header']}_RelPos_{a3_match['rel_pos']:.3f}\n{a3_match['seq']}\n")
                out.write(f">A4_{a4_match['header']}_RelPos_{a4_match['rel_pos']:.3f}\n{a4_match['seq']}\n")

print(f"Successfully clustered {group_count} high-confidence spatial trios based on true array boundaries.")
