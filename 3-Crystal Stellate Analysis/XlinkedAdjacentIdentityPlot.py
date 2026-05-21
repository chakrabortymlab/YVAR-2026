import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import os
import re

# ==========================================
# 1. BIOLOGICAL ANCHORS (VERIFIED COORDINATES)
# ==========================================
coords = [
    # ISO1
    ['ISO1', 'HetSte L1', 1, 0, 0],
    ['ISO1', 'EuSte Distal', 1, 1, 1],
    ['ISO1', 'EuSte Proximal', 10, 2, 11],
    ['ISO1', 'HetSte Singleton', 1, 12, 12],
    ['ISO1', 'HetSte L2', 4, 13, 16],
    ['ISO1', 'HetSte L3', 11, 17, 27],
    ['ISO1', 'HetSte L4', 1, 28, 28],
    ['ISO1', 'HetSte L5', 1, 29, 29],
    ['ISO1', 'HetSte L6', 1, 30, 30],
    ['ISO1', 'HetSte L7', 1, 31, 31],
    
    # A3
    ['A3', 'EuSte Distal', 1, 0, 0],
    ['A3', 'EuSte Proximal', 1, 1, 1],
    ['A3', 'HetSte L1', 1, 2, 2],
    ['A3', 'HetSte L2', 4, 3, 6],
    ['A3', 'HetSte L3', 11, 7, 17],
    ['A3', 'HetSte L4', 1, 18, 18],
    ['A3', 'HetSte L5', 1, 19, 19],
    ['A3', 'HetSte L6', 1, 20, 20],
    ['A3', 'HetSte L7', 1, 21, 21],
    
    # A4
    ['A4', 'EuSte Distal', 1, 0, 0],
    ['A4', 'EuSte Proximal Main', 188, 1, 188],
    ['A4', 'EuSte Proximal Island', 9, 189, 197],
    ['A4', 'HetSte Singleton 1', 1, 198, 198],
    ['A4', 'HetSte Singleton 2', 1, 199, 199]
]

df_coords = pd.DataFrame(coords, columns=["Strain", "Region", "Unit_Count", "Matrix_Start_Index", "Matrix_End_Index"])

# ==========================================
# 2. CONFIGURATION & AESTHETICS
# ==========================================
strains = ['ISO1', 'A3', 'A4']
gap_tolerance = 5000  

family_colors = {
    'EuSte': '#3182bd',  
    'HetSte': '#31a354', 
    'Unknown': '#d9d9d9' 
}

def parse_full_coords(header):
    match = re.search(r'([^:_#]+):(\d+)-(\d+)$', header)
    if match: return match.group(1), int(match.group(2)), int(match.group(3))
    return 'Unknown', 0, 0

# ==========================================
# 3. PLOTTING LOGIC
# ==========================================
def plot_adjacent_identities(strain):
    matrix_file = f"Stellate_Identity_{strain}_vs_{strain}.tsv"
    if not os.path.exists(matrix_file): return

    print(f"Generating clean adjacent identity plot for {strain}...")
    
    df = pd.read_csv(matrix_file, sep='\t', index_col=0)
    mat = df.values
    labels = list(df.index)
    N = len(labels)
    
    if N < 2: return
        
    adj_identities = [mat[i, i+1] for i in range(N - 1)]
    x_positions = np.arange(N - 1)
    
    fig, ax = plt.subplots(figsize=(18, 6))
    fig.patch.set_facecolor('white')
    ax.set_facecolor('#fafafa')
    
    strain_coords = df_coords[df_coords['Strain'] == strain].reset_index(drop=True)

    # Calculate dynamic Y limits to focus on variance
    min_val = min(adj_identities) if len(adj_identities) > 0 else 50
    y_min = max(50, min_val - 2)
    if y_min > 90: y_min = 90
    
    y_data_max = 100.5 
    ax.set_ylim(y_min, 114) 
    ax.set_xlim(-0.5, N - 1.5)

    # -----------------------------------------------------
    # TRACK 1: Shade Biological Domains
    # -----------------------------------------------------
    for idx, row in strain_coords.iterrows():
        start = row['Matrix_Start_Index']
        end = row['Matrix_End_Index']
        region = row['Region']
        
        btype = 'EuSte' if 'EuSte' in region else 'HetSte'
        color = family_colors.get(btype, family_colors['Unknown'])
        
        left_edge = start - 0.5
        right_edge = end + 0.5
        
        ax.fill_between([left_edge, right_edge], y_min, y_data_max, color=color, alpha=0.25, zorder=1)
        
        # Domain Boundary Line
        if right_edge < N - 1:
            ax.plot([right_edge, right_edge], [y_min, y_data_max], color='#252525', linestyle='--', linewidth=1.5, alpha=0.8, zorder=2)
            
        # Callouts
        mid_x = (left_edge + right_edge) / 2.0
        clean_label = region.replace('EuSte ', '').replace('HetSte ', '')
        
        if (end - start) > 2 or strain != 'A4':
            # Stagger text vertically to prevent overlapping 
            stagger = idx % 3
            y_line_end = y_data_max + 2 + (stagger * 3.5)
            
            # Callout Line & Anchor Dot
            ax.plot([mid_x, mid_x], [y_data_max, y_line_end], color='black', lw=1.0, zorder=5)
            ax.plot(mid_x, y_data_max, marker='o', markersize=4, color='black', zorder=6)
            
            # Label
            ax.text(mid_x, y_line_end + 0.5, clean_label, ha='center', va='bottom', rotation=0, 
                    fontsize=10, fontweight='bold', color='#333333', zorder=7,
                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=1))

    # -----------------------------------------------------
    # TRACK 2: Plot Genomic Physical Gaps (Dashed Lines Only)
    # -----------------------------------------------------
    for i in range(N - 1):
        scaf1, start1, end1 = parse_full_coords(labels[i])
        scaf2, start2, end2 = parse_full_coords(labels[i+1])
        
        if scaf1 != 'Unknown' and scaf2 != 'Unknown':
            if scaf1 != scaf2 or (start2 - end1) > gap_tolerance:
                # Draw red line only up to the data max
                ax.plot([i, i], [y_min, y_data_max], color='#cb181d', linestyle=':', linewidth=2.5, alpha=0.8, zorder=3)
    
    # -----------------------------------------------------
    # TRACK 3: Plot the Adjacent Sequence Identity Line
    # -----------------------------------------------------
    ax.plot(x_positions, adj_identities, color='#252525', linewidth=2, 
            marker='o', markersize=6, markerfacecolor='#252525', 
            markeredgecolor='white', markeredgewidth=1.2, zorder=4)

    # -----------------------------------------------------
    # STYLING & CLEANUP
    # -----------------------------------------------------
    ax.set_title(f"Sequence Identity Between Adjacent Units ({strain} X-Linked Stellate Array)", 
                 fontsize=18, fontweight='bold', pad=20, color='#111111')
    ax.set_ylabel('Sequence Identity to Next Unit (%)', fontsize=14, fontweight='bold', color='#333333')
    ax.set_xlabel('Tandem Array Position (Unit Index Transition)', fontsize=14, fontweight='bold', color='#333333')
    
    # Restrict Y-axis ticks to the data area (prevents gridlines from drawing into the header)
    tick_start = int(np.ceil(y_min / 10.0)) * 10
    ax.set_yticks(np.arange(tick_start, 101, 10))
    ax.grid(axis='y', linestyle='-', color='white', linewidth=1.5, zorder=0)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # Hide the left spine above 100%
    ax.spines['left'].set_bounds(y_min, 100)
    
    # Custom Legend securely anchored to the middle-right
    legend_elements = [
        mpatches.Patch(color='none', label='Stellate Unit Type:'),
        mpatches.Patch(facecolor=family_colors['EuSte'], alpha=0.35, edgecolor='#252525', linewidth=1, label='EuSte Array Region'),
        mpatches.Patch(facecolor=family_colors['HetSte'], alpha=0.35, edgecolor='#252525', linewidth=1, label='HetSte Array Region'),
        mpatches.Patch(color='none', label=''),
        mpatches.Patch(color='none', label='Genomic Topology:'),
        Line2D([0], [0], color='#cb181d', linestyle=':', linewidth=2.5, label='Physical Genomic Gap'),
        Line2D([0], [0], color='#252525', linestyle='--', linewidth=1.5, label='Domain Boundary')
    ]
    
    ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1.02, 0.5), 
              fontsize=12, frameon=True, facecolor='white', edgecolor='#dddddd')

    plt.tight_layout()
    plt.savefig(f"Stellate_Adjacent_Identity_{strain}.png", dpi=600, bbox_inches='tight')
    plt.savefig(f"Stellate_Adjacent_Identity_{strain}.svg", format='svg', bbox_inches='tight')
    plt.close()

# ==========================================
# EXECUTION
# ==========================================
if __name__ == "__main__":
    for strain in strains:
        plot_adjacent_identities(strain)
        
    print("\nSuccess! Final adjusted plots generated.")
