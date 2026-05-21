import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os

# ==========================================
# CONFIGURATION
# ==========================================
strains = ['ISO-1', 'A3', 'A4']

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


suste_island_sizes = {
    'ISO-1': [40, 67, 46, 27, 52, 60, 14, 52],
    'A3': [36, 39, 48, 72, 64, 48, 15, 48],
    'A4': [39, 98, 46, 73, 53, 70, 14, 50]
}

promoter_color = 'hotpink' # Distinct pink for BNACtes1

# ==========================================
# ANNOTATION LOGIC
# ==========================================
def annotate_structural_domains(lbls, strain):
    types = []
    suste_count = 0
    cumulative_thresholds = []

    if strain in suste_island_sizes:
        sizes = suste_island_sizes[strain]
        cumulative_thresholds = [sum(sizes[:i+1]) for i in range(len(sizes))]
    else:
        cumulative_thresholds = [9999] * 8 

    for lbl in lbls:
        is_pckr = 'PCKR' in lbl.upper()

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
# PLOTTING LOGIC
# ==========================================
def plot_adjacent_identities(strain, promoter_units):
    matrix_file = f"Matrix_Intra_{strain}.tsv"
    
    if not os.path.exists(matrix_file):
        print(f"Skipping {strain}: {matrix_file} not found.")
        return

    print(f"Generating adjacent identity plot for {strain}...")
    
    # Load Matrix
    df = pd.read_csv(matrix_file, sep='\t', index_col=0)
    mat = df.values
    labels = list(df.index)
    N = len(labels)
    
    # Get annotations
    types = annotate_structural_domains(labels, strain)
    
    # Extract adjacent identities (diagonal + 1)
    adj_identities = []
    for i in range(N - 1):
        adj_identities.append(mat[i, i+1])
        
    x_positions = np.arange(N - 1)
    
    # Calculate structural blocks for shading
    blocks = []
    curr_val, start = None, 0
    for i, val in enumerate(types):
        if curr_val is None: 
            curr_val = val
        elif curr_val != val:
            blocks.append((curr_val, start, i))
            curr_val = val
            start = i
    if types: 
        blocks.append((curr_val, start, len(types)))

    # Set up plot
    fig, ax = plt.subplots(figsize=(18, 6))
    fig.patch.set_facecolor('white')
    
    # Plot background shades for structural domains
    for block_type, start_idx, end_idx in blocks:
        color = array_colors.get(block_type, 'black')
        hatch = '////' if 'Inverted' in block_type else ''
        alpha_val = 0.3 if 'PCKR' not in block_type else 0.15
        
        ax.axvspan(start_idx, end_idx, color=color, alpha=alpha_val, hatch=hatch, lw=0)
        
        # Draw hard vertical lines at boundaries
        if end_idx < N:
            ax.axvline(x=end_idx, color='black', linestyle='--', linewidth=1.5, alpha=0.7)

    # Plot BNACtes1 Promoters completely ABOVE the plot axis
    for i, lbl in enumerate(labels):
        if lbl in promoter_units:
            # Using axes transform for Y (1.0 = top line of plot, so 1.02 to 1.06 hovers above it)
            # Using data coordinates for X (i = the unit index)
            # clip_on=False ensures it is allowed to draw outside the graph borders
            ax.plot([i, i], [1.02, 1.03], transform=ax.get_xaxis_transform(), 
                    color=promoter_color, linestyle='-', linewidth=2, alpha=1.0, clip_on=False, zorder=5)

    # Plot the adjacent identity line (on top of the blocks)
    ax.plot(x_positions, adj_identities, color='black', linewidth=1.5, marker='o', markersize=3, markerfacecolor='red', markeredgewidth=0, zorder=4)

    # Styling
    ax.set_title(f"Sequence Identity Drops at Structural Boundaries ({strain} Array)", fontsize=18, fontweight='bold', pad=30) 
    ax.set_ylabel('Sequence Identity to Next Unit (%)', fontsize=14, fontweight='bold')
    ax.set_xlabel('Tandem Array Position (Unit Index)', fontsize=14, fontweight='bold')
    
    # Set limits to zoom in on the critical range
    min_val = min(adj_identities)
    ax.set_ylim(max(50, min_val - 2), 100.5) 
    ax.set_xlim(0, N-1)
    
    ax.grid(axis='y', linestyle=':', alpha=0.6)
    
    # Add a custom legend
    legend_elements = [
        mpatches.Patch(color='none', label='Structural Domains:'),
        mpatches.Patch(facecolor=array_colors['PCKR_Forward'], alpha=0.3, edgecolor='black', label='PCKR (Forward)'),
        mpatches.Patch(facecolor=array_colors['PCKR_Inverted'], alpha=0.3, edgecolor='black', hatch='////', label='PCKR (Inverted Anchors)')
    ]
    for i in range(1, 8):
        legend_elements.append(mpatches.Patch(facecolor=array_colors[f'SuSte_{i}'], alpha=0.3, edgecolor='black', label=f'SuSte (Island {i})'))
    legend_elements.append(mpatches.Patch(facecolor=array_colors['SuSte_8_Inverted'], alpha=0.3, edgecolor='black', hatch='////', label='SuSte (Island 8, Inverted)'))
    
    # Add the Promoter to the legend
    legend_elements.append(mpatches.Patch(color='none', label=''))
    legend_elements.append(mpatches.Patch(color='none', label='Features:'))
    legend_elements.append(mpatches.Patch(facecolor=promoter_color, alpha=1.0, edgecolor='none', label='BNACtes1 Promoter (Above Plot)'))

    # Place legend outside to the right
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.01, 1), fontsize=10, frameon=False)

    plt.tight_layout()
    plt.savefig(f"Adjacent_Identity_Drops_{strain}.png", dpi=600, bbox_inches='tight')
    plt.savefig(f"Adjacent_Identity_Drops_{strain}.svg", format='svg', dpi=600, bbox_inches='tight')
    plt.close()

# ==========================================
# EXECUTION
# ==========================================
# 1. Load Promoter Annotations
promoter_units = set()
bnactes_file = 'BNACtes_Unit_Locations.tsv'
if os.path.exists(bnactes_file):
    print(f"Loading BNACtes Annotations from {bnactes_file}...")
    df_prom = pd.read_csv(bnactes_file, sep='\t')
    promoter_units = set(df_prom['Unit_Header'])
else:
    print(f"Warning: {bnactes_file} not found. Promoter tracks will be empty.")

# 2. Generate Plots
for strain in strains:
    plot_adjacent_identities(strain, promoter_units)
    
print("Success! Adjacent Sequence Identity plots generated.")
