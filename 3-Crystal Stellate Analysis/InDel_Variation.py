import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patheffects as path_effects

# ==========================================
# 1. FASTA PARSING & GLOBAL CONSENSUS BUILDER
# ==========================================
def parse_fasta(file_path):
    headers = []
    seqs = []
    with open(file_path, 'r') as f:
        curr_seq = []
        for line in f:
            if line.startswith('>'):
                headers.append(line.strip())
                if curr_seq:
                    seqs.append(''.join(curr_seq))
                    curr_seq = []
            else:
                curr_seq.append(line.strip().upper())
        if curr_seq:
            seqs.append(''.join(curr_seq))
    return headers, seqs

def process_alignment(fasta_path, strain, window=10, phase_start=None):
    headers, all_seqs = parse_fasta(fasta_path)
    
    if not all_seqs:
        return pd.DataFrame()
        
    # --- 1. BUILD THE GLOBAL CONSENSUS BACKBONE ---
    max_len = max(len(s) for s in all_seqs)
    padded_all = [s.ljust(max_len, '-') for s in all_seqs]
    arr_all = np.array([list(s) for s in padded_all])
    
    total_seqs, align_len = arr_all.shape
    global_gap_pcts = np.sum(arr_all == '-', axis=0) / total_seqs * 100
    is_cons_col = global_gap_pcts < 50 # The Global Mask
    
    # --- 2. FILTER FOR THE SPECIFIC STRAIN ---
    strain_indices = [i for i, h in enumerate(headers) if strain in h]
    if not strain_indices:
        print(f"  -> No sequences found for strain {strain} in {fasta_path}. Skipping.")
        return pd.DataFrame()
        
    arr_strain = arr_all[strain_indices, :]
    num_seqs = arr_strain.shape[0]
    
    current_cons_pos = 0
    cons_data = {}
    seqs_with_insertion = np.zeros(num_seqs, dtype=bool)
    
    # --- 3. CALCULATE METRICS USING GLOBAL BACKBONE ---
    for col_idx in range(align_len):
        col = arr_strain[:, col_idx]
        is_cons = is_cons_col[col_idx] # Crucial: using the global boolean mask
        
        if is_cons:
            current_cons_pos += 1
            
            gaps = np.sum(col == '-')
            del_pct = (gaps / num_seqs) * 100
            
            chars, counts = np.unique(col, return_counts=True)
            valid_mask = chars != '-'
            valid_chars = chars[valid_mask]
            valid_counts = counts[valid_mask]
            
            if len(valid_chars) > 0:
                cons_char = valid_chars[np.argmax(valid_counts)]
                snps = np.sum((col != '-') & (col != cons_char))
                snp_pct = (snps / num_seqs) * 100
            else:
                snp_pct = 0.0
                
            cons_data[current_cons_pos] = {
                'Consensus_Pos': current_cons_pos,
                'Deletion_Pct': del_pct,
                'SNP_Pct': snp_pct,
                'Insertion_Pct': 0.0
            }
            seqs_with_insertion = np.zeros(num_seqs, dtype=bool)
            
        else:
            has_insert = col != '-'
            seqs_with_insertion = seqs_with_insertion | has_insert
            
            pos_to_assign = max(1, current_cons_pos)
            if pos_to_assign in cons_data:
                ins_pct = (np.sum(seqs_with_insertion) / num_seqs) * 100
                cons_data[pos_to_assign]['Insertion_Pct'] = ins_pct
                
    df = pd.DataFrame(list(cons_data.values()))
    
    # --- 4. PHASE SHIFTING LOGIC ---
    if not df.empty and phase_start is not None:
        L = df['Consensus_Pos'].max()
        shift_amount = phase_start - 1
        
        df['Consensus_Pos'] = df['Consensus_Pos'] - shift_amount
        df.loc[df['Consensus_Pos'] <= 0, 'Consensus_Pos'] += L
        
        df = df.sort_values('Consensus_Pos').reset_index(drop=True)
    
    if not df.empty:
        df['Del_Smooth'] = df['Deletion_Pct'].rolling(window, center=True, min_periods=1).mean()
        df['SNP_Smooth'] = df['SNP_Pct'].rolling(window, center=True, min_periods=1).mean()
        df['Ins_Smooth'] = df['Insertion_Pct'].rolling(window, center=True, min_periods=1).mean()
    return df

# ==========================================
# 2. ANNOTATION SHIFTING HELPER
# ==========================================
def shift_annotations(anns, phase_start, total_len):
    """Circularly shifts annotation blocks while preserving custom colors/attributes."""
    new_anns = []
    shift_amount = phase_start - 1
    
    for ann in anns:
        s = ann['start'] - shift_amount
        e = ann['end'] - shift_amount
        
        if s <= 0: s += total_len
        if e <= 0: e += total_len
        
        ann_p1 = ann.copy()
        ann_p2 = ann.copy()
        
        if s > e:
            ann_p1.update({'start': s, 'end': total_len})
            ann_p2.update({'start': 1, 'end': e})
            new_anns.extend([ann_p1, ann_p2])
        else:
            ann_p1.update({'start': s, 'end': e})
            new_anns.append(ann_p1)
            
    return new_anns

# ==========================================
# 3. PLOTTING & GENE TRACK ANNOTATION
# ==========================================
def plot_clean_profile(ax, df, color_del, color_snp, color_ins, title, annotations):
    if df.empty: return
    x = df['Consensus_Pos']
    y_del = df['Del_Smooth']
    y_snp = df['SNP_Smooth']
    y_ins = df['Ins_Smooth']
    
    max_y = df[['Del_Smooth', 'SNP_Smooth', 'Ins_Smooth']].max().max()
    if pd.isna(max_y) or max_y < 5: 
        max_y = 5 
        
    y_rail = -0.10 * max_y          
    y_rect_bottom = -0.14 * max_y   
    rect_height = 0.08 * max_y      
    y_text = -0.18 * max_y          
    
    ax.set_ylim(-0.25 * max_y, max_y * 1.1)
    ax.set_xlim(0, x.max())
    
    ax.fill_between(x, 0, y_del, color=color_del, alpha=0.3)
    ax.plot(x, y_del, color=color_del, linewidth=2.5, label='% Deletion (Missing Sequence)')
    
    ax.fill_between(x, 0, y_snp, color=color_snp, alpha=0.3)
    ax.plot(x, y_snp, color=color_snp, linewidth=2.5, linestyle='--', label='% Divergence (SNPs/Mismatches)')
    
    ax.fill_between(x, 0, y_ins, color=color_ins, alpha=0.2)
    ax.plot(x, y_ins, color=color_ins, linewidth=2.5, linestyle=':', label='% Insertion (Extra Sequence)')
    
    ax.axhline(0, color='black', linewidth=1.5, zorder=3)
    ax.plot([0, x.max()], [y_rail, y_rail], color='#e2e8f0', linewidth=6, zorder=1)
    
    colors = ['#cbd5e1', '#94a3b8', '#64748b', '#475569']
    for i, ann in enumerate(annotations):
        start = ann['start']
        end = ann['end']
        label = ann['label']
        c = ann.get('color', colors[i % len(colors)]) 
        
        rect = Rectangle((start, y_rect_bottom), end - start, rect_height, facecolor=c, edgecolor='black', linewidth=1, zorder=2)
        ax.add_patch(rect)
        
        txt = ax.text((start + end)/2, y_text, label, ha='center', va='top', fontsize=10, fontweight='bold', color='#1e293b')
        txt.set_path_effects([path_effects.withStroke(linewidth=2, foreground='white')])
        
        ax.axvspan(start, end, color='gray', alpha=0.08, zorder=0)

    ax.set_title(title, fontsize=18, fontweight='bold', pad=15)
    ax.set_ylabel('% of Copies in Array', fontsize=14, fontweight='bold')
    
    ticks = ax.get_yticks()
    ax.set_yticks([t for t in ticks if t >= 0])
    
    ax.grid(axis='y', linestyle='--', alpha=0.5)
    ax.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.9, fontsize=11)

# Define the regions (with Satellites included)
pckr_anns = [
    {'start': 308, 'end': 1631, 'label': 'HeT-A Relic'},
    {'start': 1777, 'end': 1820, 'label': 'A-rich', 'color': '#abd9e9'},
    {'start': 1867, 'end': 1985, 'label': 'Copia LTR'},
    {'start': 2263, 'end': 2563, 'label': 'CG40635'}
]

suste_anns = [
    {'start': 50, 'end': 1105, 'label': 'Hoppel / PROTOP\nInsertion'},
    {'start': 1099, 'end': 1149, 'label': 'βNACtes1'},
    {'start': 1472, 'end': 1694, 'label': 'A-rich / Simple', 'color': '#abd9e9'},
    {'start': 2166, 'end': 2651, 'label': 'Stellate\nHomology'},
    {'start': 2707, 'end': 2763, 'label': 'βNACtes1'}
]

# ==========================================
# 4. GENERATE THE FIGURES (SVG & PNG)
# ==========================================
strains = ['ISO-1', 'A3', 'A4']

print("Processing MAFFT alignments and generating SVG/PNG files...")

for strain in strains:
    print(f"\n--- Processing Strain: {strain} ---")
    
    # ------------------
    # PCKR Figures
    # ------------------
    df_pckr = process_alignment('PCKR_AlignedMAFFT.fasta', strain, window=10)
    
    if not df_pckr.empty:
        fig, ax = plt.subplots(figsize=(16, 6))
        fig.patch.set_facecolor('white')
        
        plot_clean_profile(
            ax, df_pckr, '#e11d48', '#f59e0b', '#9333ea', 
            f'PCKR Array ({strain}): Structural & Sequence Variation Profile', 
            pckr_anns
        )
        ax.set_xlabel('Global Consensus Coordinate (bp)', fontsize=16, fontweight='bold')
        plt.tight_layout()
        
        output_base = f'PCKR_{strain}_Variation_Profile'
        plt.savefig(f'{output_base}.png', dpi=300, bbox_inches='tight')
        plt.savefig(f'{output_base}.svg', format='svg', bbox_inches='tight')
        plt.close()
        print(f"Saved {output_base}.png and .svg")

    # ------------------
    # SuSte Figures
    # ------------------
    phase_start_suste = 2707 
    df_suste = process_alignment('SuSte_AlignedMAFFT.fasta', strain, window=10, phase_start=phase_start_suste)
    
    if not df_suste.empty:
        total_len = df_suste['Consensus_Pos'].max()
        phased_suste_anns = shift_annotations(suste_anns, phase_start_suste, total_len)
        
        fig, ax = plt.subplots(figsize=(16, 6))
        fig.patch.set_facecolor('white')
        
        plot_clean_profile(
            ax, df_suste, '#0284c7', '#10b981', '#8b5cf6', 
            f'Su(Ste) Array ({strain}): Structural & Sequence Variation Profile (Phased)', 
            phased_suste_anns
        )
        ax.set_xlabel('Global Consensus Coordinate (bp)', fontsize=16, fontweight='bold')
        plt.tight_layout()
        
        output_base = f'SuSte_{strain}_Variation_Profile'
        plt.savefig(f'{output_base}.png', dpi=300, bbox_inches='tight')
        plt.savefig(f'{output_base}.svg', format='svg', bbox_inches='tight')
        plt.close()
        print(f"Saved {output_base}.png and .svg")

print("\nAll figures processed successfully!")
