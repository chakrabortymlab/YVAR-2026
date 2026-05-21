import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd
import os
import sys
import re

# ==========================================
# CONFIGURATION
# ==========================================
aln_file = 'PCKR_SuSte_Combined_Aligned.fasta'
strains = ['ISO-1', 'A3', 'A4']

min_identity_scale = 50 

okabe_ito_gradient = [
    "#0072B2", # Dark Blue (Bottom end)
    "#56B4E9", # Sky Blue
    "#009E73", # Bluish Green
    "#F0E442", # Yellow (starts around the 85-90% mark)
    "#E69F00", # Orange (Mid 90s)
    "#8C3E00"  # Darkened Vermilion / Deep Rust (100%) - Provides sharp contrast against the yellow/orange!
]
okabe_ito_cmap = LinearSegmentedColormap.from_list("OkabeIto_HighContrast", okabe_ito_gradient)

array_colors = {
    'PCKR_Forward': '#2ca02c',       
    'PCKR_Inverted': '#2ca02c',      
    'SuSte_1': '#c6dbef',            
    'SuSte_2': '#9ecae1',            
    'SuSte_3': '#6baed6',            
    'SuSte_4': '#4292c6',            
    'SuSte_5': '#2171b5',            
    'SuSte_6': '#08519c',            
    'SuSte_7': '#08306b',            
    'SuSte_8_Inverted': '#041836'    
}

strain_colors = {'ISO-1': '#9467bd', 'A3': '#ff7f0e', 'A4': '#17becf'} 
promoter_color = '#e377c2' # Tab:Pink for BNACtes1

suste_island_sizes = {
    'ISO-1': [40, 67, 46, 27, 52, 60, 14, 52],
    'A3': [36, 39, 48, 72, 64, 48, 15, 48],
    'A4': [39, 98, 46, 73, 53, 70, 14, 50]
}

# ==========================================
# 1. PARSE & SORT FASTA 
# ==========================================
def read_and_sort_fasta(filename):
    if not os.path.exists(filename):
        print(f"Error: {filename} not found.")
        sys.exit()
        
    sequences = {}
    current_name, current_seq = None, []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_name: sequences[current_name] = "".join(current_seq)
                current_name = line[1:]
                current_seq = []
            else:
                current_seq.append(line.upper())
    if current_name: sequences[current_name] = "".join(current_seq)

    def get_sort_key(header):
        match = re.match(r'^(.*)_([^:]+):(\d+)-(\d+)$', header)
        if match:
            prefix, scaf, start, end = match.groups()
            return (scaf, int(start))
        return ('Unknown', 0)

    strain_seqs = {s: {} for s in strains}
    for header, seq in sequences.items():
        for s in strains:
            if header.startswith(s):
                strain_seqs[s][header] = seq
                break
                
    sorted_strain_seqs = {}
    for s in strains:
        sorted_keys = sorted(strain_seqs[s].keys(), key=get_sort_key)
        sorted_strain_seqs[s] = {k: strain_seqs[s][k] for k in sorted_keys}
        
    return sorted_strain_seqs

# ==========================================
# 2. MATRIX MATH ENGINE 
# ==========================================
def calculate_identity_matrix(dict_seqs):
    names = list(dict_seqs.keys())
    arr = [np.frombuffer(s.encode('ascii'), dtype=np.uint8) for s in dict_seqs.values()]
    N = len(names)
    matrix = np.zeros((N, N))
    hyphen = ord('-')
    
    for i in range(N):
        for j in range(i, N): 
            sy = arr[i]
            sx = arr[j]
            either_gap = (sy == hyphen) | (sx == hyphen)
            match = (sy == sx) & ~either_gap
            matches = np.sum(match)
            valid_len = len(sy) - np.sum(either_gap) 
            
            val = (matches / valid_len * 100) if valid_len > 0 else 0
            matrix[i, j] = val
            matrix[j, i] = val 
            
    return matrix, names

# ==========================================
# 3. STRUCTURAL ANNOTATION LOGIC
# ==========================================
def annotate_structural_domains(lbls):
    types = []
    current_strain = None
    suste_count = 0
    cumulative_thresholds = []

    for lbl in lbls:
        strain = lbl.split('_')[0]
        is_pckr = 'PCKR' in lbl.upper()

        if strain != current_strain:
            current_strain = strain
            suste_count = 0
            if strain in suste_island_sizes:
                sizes = suste_island_sizes[strain]
                cumulative_thresholds = [sum(sizes[:i+1]) for i in range(len(sizes))]
            else:
                cumulative_thresholds = [9999] * 8 

        if is_pckr:
            if suste_count == 0:
                types.append('PCKR_Forward')
            else:
                types.append('PCKR_Inverted') 
        else:
            suste_count += 1
            island_idx = 1
            for i, thresh in enumerate(cumulative_thresholds):
                if suste_count <= thresh:
                    island_idx = i + 1
                    break
            
            if island_idx >= 8:
                types.append('SuSte_8_Inverted') 
            else:
                types.append(f'SuSte_{island_idx}')

    return types

# ==========================================
# 4. EXPORT ISLAND TSV DATA
# ==========================================
def export_island_data(strain, matrix, labels, types):
    print(f"Exporting Structural Data TSVs for {strain}...")
    summary_data = []
    unit_data = []
    
    unique_types = []
    for t in types:
        if t not in unique_types:
            unique_types.append(t)
            
    def get_coords(h):
        match = re.search(r':(\d+)-(\d+)$', h)
        return int(match.group(1)), int(match.group(2)) if match else (0, 0)
            
    for domain in unique_types:
        indices = [i for i, x in enumerate(types) if x == domain]
        size = len(indices)
        if size == 0: continue
        
        first_header = labels[indices[0]]
        last_header = labels[indices[-1]]
        
        start_coord, _ = get_coords(first_header)
        _, end_coord = get_coords(last_header)
        
        # Calculate intra-domain mean identity
        if size > 1:
            sub_mat = matrix[np.ix_(indices, indices)]
            intra_id = np.mean(sub_mat[np.triu_indices(size, k=1)])
        else:
            intra_id = 100.0
            
        summary_data.append({
            'Strain': strain,
            'Structural_Domain': domain,
            'Copy_Number': size,
            'Start_Coordinate': start_coord,
            'End_Coordinate': end_coord,
            'Intra_Domain_Identity_Pct': round(intra_id, 3)
        })
        
        # Detailed Unit Stats
        for idx in indices:
            header = labels[idx]
            s_coord, e_coord = get_coords(header)
            
            if size > 1:
                row_vals = matrix[idx, indices]
                row_vals = np.delete(row_vals, indices.index(idx)) # remove self-comparison
                unit_mean = np.mean(row_vals)
            else:
                unit_mean = 100.0
                
            unit_data.append({
                'Strain': strain,
                'Structural_Domain': domain,
                'Unit_Header': header,
                'Start_Coordinate': s_coord,
                'End_Coordinate': e_coord,
                'Mean_Identity_To_Domain_Pct': round(unit_mean, 3)
            })
            
    df_summary = pd.DataFrame(summary_data)
    df_units = pd.DataFrame(unit_data)
    df_summary.to_csv(f"Island_Summary_{strain}.tsv", sep='\t', index=False)
    df_units.to_csv(f"Detailed_Unit_Stats_{strain}.tsv", sep='\t', index=False)

# ==========================================
# 5. HEATMAP PLOTTING
# ==========================================
def plot_triangle_heatmap(matrix, labels, types, title, filename, promoter_units=None):
    if promoter_units is None: promoter_units = set()
    print(f"Plotting {filename}...")
    N = len(labels)
    fig, ax = plt.subplots(figsize=(16, 10))
    fig.patch.set_facecolor('white')

    i_grid, j_grid = np.meshgrid(np.arange(N+1), np.arange(N+1), indexing='ij')
    X = (i_grid + j_grid) / 2.0
    Y = (j_grid - i_grid) / 2.0

    C = np.where(i_grid[:-1, :-1] <= j_grid[:-1, :-1], matrix, np.nan)
    
    mesh = ax.pcolormesh(X, Y, C, cmap=okabe_ito_cmap, vmin=min_identity_scale, vmax=100, rasterized=True)

    cbar = plt.colorbar(mesh, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label('Sequence Identity (%)', rotation=270, labelpad=25, fontsize=14, fontweight='bold')
    cbar.ax.tick_params(labelsize=12)
    
    ax.set_aspect('equal')

    # Visual Blocking setup
    array_blocks = []
    curr_val, start = None, 0
    for i, val in enumerate(types):
        if curr_val is None: curr_val = val
        elif curr_val != val:
            array_blocks.append((curr_val, start, i))
            curr_val = val
            start = i
    if types: array_blocks.append((curr_val, start, len(types)))

    strain_blocks = []
    curr_strain, start = None, 0
    for i, lbl in enumerate(labels):
        val = lbl.split('_')[0]
        if curr_strain is None: curr_strain = val
        elif curr_strain != val:
            strain_blocks.append((curr_strain, start, i))
            curr_strain = val
            start = i
    if labels: strain_blocks.append((curr_strain, start, len(labels)))

    for strain, start, end in strain_blocks:
        if start > 0:
            b = start
            ax.plot([b, (b + N)/2], [0, (N - b)/2], color='white', linestyle='--', linewidth=2)
            ax.plot([b/2, b], [b/2, 0], color='white', linestyle='--', linewidth=2)

    # --- ENLARGED TRACK COEFFS (TRIANGLE) ---
    h_prom = N * 0.035    # Increased from 0.010
    y_prom = -N * 0.040   # Offset rearranged to clear bottom matrix baseline
    
    h_array = N * 0.045   # Increased from 0.020
    y_array = -N * 0.090  # Offset adjusted down
    
    h_strain = N * 0.030  # Increased from 0.015
    y_strain = -N * 0.130 # Offset adjusted down

    # 1. Promoter Track (Topmost track near the matrix)
    for i, lbl in enumerate(labels):
        if lbl in promoter_units:
            rect = mpatches.Rectangle((i, y_prom), width=1, height=h_prom, 
                                      facecolor=promoter_color, edgecolor='none', lw=0)
            ax.add_patch(rect)

    # 2. Structural Array Track
    for block_type, start, end in array_blocks:
        color = array_colors.get(block_type, 'black')
        hatch_pattern = '////' if 'Inverted' in block_type else ''
        rect = mpatches.Rectangle((start, y_array), width=(end-start), height=h_array, 
                                  facecolor=color, edgecolor='black', lw=0.8, hatch=hatch_pattern)
        ax.add_patch(rect)

    # 3. Strain Track
    if len(strain_blocks) > 1:
        for strain, start, end in strain_blocks:
            s_color = strain_colors.get(strain, 'black')
            rect = mpatches.Rectangle((start, y_strain), width=(end-start), height=h_strain, 
                                      facecolor=s_color, edgecolor='black', lw=0.8)
            ax.add_patch(rect)

    ax.set_title(title, fontsize=18, fontweight='bold', pad=20)
    ax.axis('off') 
    ax.set_ylim(bottom=-N*0.16, top=N*0.55) # Extensively expanded bottom headroom limit from -0.09 to -0.16
    ax.set_xlim(left=-N*0.05, right=N*1.05)

    plt.tight_layout()
    plt.savefig(f"{filename}.png", dpi=600)
    plt.savefig(f"{filename}.svg", format='svg', dpi=600)
    plt.close()

def plot_rectangular_heatmap(sub_matrix, types_y, types_x, strain_y, strain_x, labels_y, labels_x, title, filename, promoter_units=None):
    if promoter_units is None: promoter_units = set()
    print(f"Plotting {filename}...")
    Ny, Nx = sub_matrix.shape
    fig, ax = plt.subplots(figsize=(16, 10))
    fig.patch.set_facecolor('white')

    X, Y = np.meshgrid(np.arange(Nx+1), np.arange(Ny+1))

    mesh = ax.pcolormesh(X, Y, sub_matrix, cmap=okabe_ito_cmap, vmin=min_identity_scale, vmax=100, rasterized=True)

    cbar = plt.colorbar(mesh, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label('Sequence Identity (%)', rotation=270, labelpad=25, fontsize=14, fontweight='bold')
    cbar.ax.tick_params(labelsize=12)

    ax.set_aspect('equal')

    # X axis blocks (Bottom)
    array_blocks_x = []
    curr_val, start = None, 0
    for i, val in enumerate(types_x):
        if curr_val is None: curr_val = val
        elif curr_val != val:
            array_blocks_x.append((curr_val, start, i))
            curr_val = val
            start = i
    if types_x: array_blocks_x.append((curr_val, start, len(types_x)))

    # Y axis blocks (Left)
    array_blocks_y = []
    curr_val, start = None, 0
    for i, val in enumerate(types_y):
        if curr_val is None: curr_val = val
        elif curr_val != val:
            array_blocks_y.append((curr_val, start, i))
            curr_val = val
            start = i
    if types_y: array_blocks_y.append((curr_val, start, len(types_y)))

    # --- ENLARGED TRACK COEFFS (RECTANGULAR X-AXIS & Y-AXIS) ---
    h_prom_x = Ny * 0.035   # Increased from 0.010
    y_prom_x = -Ny * 0.040  # Reset offset position
    h_array_x = Ny * 0.045  # Increased from 0.020
    y_pos_x = -Ny * 0.090   # Reset offset position

    w_prom_y = Nx * 0.035   # Increased from 0.010
    x_prom_y = -Nx * 0.040  # Reset offset position
    w_array_y = Nx * 0.045  # Increased from 0.020
    x_pos_y = -Nx * 0.090   # Reset offset position

    # X-Axis Promoter Track
    for i, lbl in enumerate(labels_x):
        if lbl in promoter_units:
            rect = mpatches.Rectangle((i, y_prom_x), width=1, height=h_prom_x, 
                                      facecolor=promoter_color, edgecolor='none', clip_on=False)
            ax.add_patch(rect)

    # X-Axis Structural Track
    for block_type, start, end in array_blocks_x:
        color = array_colors.get(block_type, 'black')
        hatch_pattern = '////' if 'Inverted' in block_type else ''
        rect = mpatches.Rectangle((start, y_pos_x), width=(end-start), height=h_array_x, 
                                  facecolor=color, edgecolor='black', lw=0.8, hatch=hatch_pattern, clip_on=False)
        ax.add_patch(rect)

    # Y-Axis Promoter Track
    for i, lbl in enumerate(labels_y):
        if lbl in promoter_units:
            rect = mpatches.Rectangle((x_prom_y, i), width=w_prom_y, height=1, 
                                      facecolor=promoter_color, edgecolor='none', clip_on=False)
            ax.add_patch(rect)

    # Y-Axis Structural Track
    for block_type, start, end in array_blocks_y:
        color = array_colors.get(block_type, 'black')
        hatch_pattern = '////' if 'Inverted' in block_type else ''
        rect = mpatches.Rectangle((x_pos_y, start), width=w_array_y, height=(end-start), 
                                  facecolor=color, edgecolor='black', lw=0.8, hatch=hatch_pattern, clip_on=False)
        ax.add_patch(rect)

    ax.set_title(title, fontsize=18, fontweight='bold', pad=20)
    
    # Expanded labelpad limits from 0.06 to 0.14 so text doesn't overlap tracks
    ax.set_xlabel(f"{strain_x} Array Units", fontsize=14, fontweight='bold', labelpad=Ny*0.14)
    ax.set_ylabel(f"{strain_y} Array Units", fontsize=14, fontweight='bold', labelpad=Nx*0.14)
    ax.set_xticks([])
    ax.set_yticks([])
    
    for spine in ax.spines.values():
        spine.set_visible(False)

    plt.tight_layout()
    plt.savefig(f"{filename}.png", dpi=600, bbox_inches='tight')
    plt.savefig(f"{filename}.svg", format='svg', dpi=600, bbox_inches='tight')
    plt.close()

# ==========================================
# 6. STANDALONE LEGEND
# ==========================================
def plot_standalone_legend():
    print("Generating Combined Master Legend...")
    fig, ax = plt.subplots(figsize=(6, 8))
    fig.patch.set_facecolor('white')
    ax.axis('off')

    legend_elements = []
    legend_elements.append(mpatches.Patch(color='none', label='Structural Domains:'))
    legend_elements.append(mpatches.Patch(facecolor=array_colors['PCKR_Forward'], edgecolor='black', label='  PCKR (Forward Contig)'))
    legend_elements.append(mpatches.Patch(facecolor=array_colors['PCKR_Inverted'], edgecolor='black', hatch='////', label='  PCKR (Inverted Anchors)'))
    
    for i in range(1, 8):
        legend_elements.append(mpatches.Patch(facecolor=array_colors[f'SuSte_{i}'], edgecolor='black', label=f'  SuSte (Island {i})'))
    legend_elements.append(mpatches.Patch(facecolor=array_colors['SuSte_8_Inverted'], edgecolor='black', hatch='////', label='  SuSte (Island 8, Inverted)'))
        
    legend_elements.append(mpatches.Patch(color='none', label='')) 
    legend_elements.append(mpatches.Patch(color='none', label='Features:'))
    legend_elements.append(mpatches.Patch(facecolor=promoter_color, edgecolor='none', label='  BNACtes1 Promoter'))

    legend_elements.append(mpatches.Patch(color='none', label='')) 
    legend_elements.append(mpatches.Patch(color='none', label='Strain Map:'))
    for k, v in strain_colors.items():
        legend_elements.append(mpatches.Patch(facecolor=v, edgecolor='black', label=f'  {k}'))

    ax.legend(handles=legend_elements, loc='center', fontsize=12, frameon=False, title_fontsize=14, labelspacing=0.8)
    plt.tight_layout()
    plt.savefig("Combined_Master_Legend.png", dpi=600)
    plt.savefig("Combined_Master_Legend.svg", format='svg', dpi=600)
    plt.close()

# ==========================================
# EXECUTION
# ==========================================
# Load Promoter Annotations
promoter_units = set()
bnactes_file = 'BNACtes_Unit_Locations_Updated.tsv'
if os.path.exists(bnactes_file):
    print(f"Loading BNACtes Annotations from {bnactes_file}...")
    df_prom = pd.read_csv(bnactes_file, sep='\t')
    promoter_units = set(df_prom['Unit_Header'])
else:
    print(f"Warning: {bnactes_file} not found. Promoter track will be empty.")

print("Loading aligned sequences...")
seqs = read_and_sort_fasta(aln_file)

# 1. INTRA-STRAIN
for strain in strains:
    print(f"\nProcessing Intra-Strain: {strain}")
    mat, lbls = calculate_identity_matrix(seqs[strain])
    types = annotate_structural_domains(lbls)
    
    export_island_data(strain, mat, lbls, types)
    
    pd.DataFrame(mat, index=lbls, columns=lbls).to_csv(f"Matrix_Intra_{strain}.tsv", sep='\t')
    plot_triangle_heatmap(mat, lbls, types, f"{strain} Contiguous Array Divergence", f"Triangle_Intra_{strain}", promoter_units=promoter_units)

# 2. INTER-STRAIN
inter_pairs = [('A4', 'ISO-1'), ('A4', 'A3'), ('ISO-1', 'A3')]
for s1, s2 in inter_pairs:
    print(f"\nProcessing Inter-Strain: {s1} vs {s2}")
    combined_seqs = {**seqs[s1], **seqs[s2]} 
    mat, lbls = calculate_identity_matrix(combined_seqs)
    types = annotate_structural_domains(lbls)
    
    pd.DataFrame(mat, index=lbls, columns=lbls).to_csv(f"Matrix_Inter_{s1}_vs_{s2}.tsv", sep='\t')
    plot_triangle_heatmap(mat, lbls, types, f"Evolutionary Divergence: {s1} vs {s2}", f"Triangle_Inter_{s1}_vs_{s2}", promoter_units=promoter_units)

    # Rectangular Heatmap Execution
    len_s1 = len(seqs[s1])
    
    # Extract the top-right block matrix (Strain 1 vs Strain 2 exclusively)
    sub_mat = mat[:len_s1, len_s1:]
    
    types_s1 = types[:len_s1]
    types_s2 = types[len_s1:]
    labels_s1 = lbls[:len_s1]
    labels_s2 = lbls[len_s1:]
    
    plot_rectangular_heatmap(
        sub_matrix=sub_mat, 
        types_y=types_s1, 
        types_x=types_s2, 
        strain_y=s1, 
        strain_x=s2, 
        labels_y=labels_s1,
        labels_x=labels_s2,
        title=f"Orthologous Array Comparison: {s1} vs {s2}", 
        filename=f"Rectangular_Inter_{s1}_vs_{s2}",
        promoter_units=promoter_units
    )

# 3. LEGEND
plot_standalone_legend()

print("\nSuccess! All SVGs, PNGs, TSVs, and Annotated Data Files have been generated.")
