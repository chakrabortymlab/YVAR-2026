import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import os

# ==========================================
# 1. BIOLOGICAL ANCHORS (VERIFIED COORDINATES)
# ==========================================
coords = [
    # ISO1 (32 total units)
    ['ISO1', 'HetSte L1', 1, 0, 0],
    ['ISO1', 'EuSte Distal', 1, 1, 1],
    ['ISO1', 'EuSte Proximal', 10, 2, 11],          # Corrected to end at Index 11
    ['ISO1', 'HetSte Singleton', 1, 12, 12],        # Separated Ste_13_HET into its own block
    ['ISO1', 'HetSte L2', 4, 13, 16],
    ['ISO1', 'HetSte L3', 11, 17, 27],
    ['ISO1', 'HetSte L4', 1, 28, 28],
    ['ISO1', 'HetSte L5', 1, 29, 29],
    ['ISO1', 'HetSte L6', 1, 30, 30],
    ['ISO1', 'HetSte L7', 1, 31, 31],
    
    # A3 (22 total units)
    ['A3', 'EuSte Distal', 1, 0, 0],
    ['A3', 'EuSte Proximal', 1, 1, 1],
    ['A3', 'HetSte L1', 1, 2, 2],
    ['A3', 'HetSte L2', 4, 3, 6],
    ['A3', 'HetSte L3', 11, 7, 17],
    ['A3', 'HetSte L4', 1, 18, 18],
    ['A3', 'HetSte L5', 1, 19, 19],
    ['A3', 'HetSte L6', 1, 20, 20],
    ['A3', 'HetSte L7', 1, 21, 21],
    
    # A4 (200 total units)
    ['A4', 'EuSte Distal', 1, 0, 0],
    ['A4', 'EuSte Proximal Main', 188, 1, 188],
    ['A4', 'EuSte Proximal Island', 9, 189, 197],
    ['A4', 'HetSte Singleton 1', 1, 198, 198],
    ['A4', 'HetSte Singleton 2', 1, 199, 199]
]

df_coords = pd.DataFrame(coords, columns=["Strain", "Region", "Unit_Count", "Matrix_Start_Index", "Matrix_End_Index"])
df_coords.to_csv("Heatmap_Region_Coordinates_Verified.tsv", sep='\t', index=False)

# ==========================================
# 2. CONFIGURATION & STYLING
# ==========================================
strains = ['ISO1', 'A3', 'A4']
min_identity_scale = 60

family_colors = {
    'EuSte': '#3182bd',  
    'HetSte': '#31a354'  
}

def add_figure_legend(ax):
    legend_elements = [
        mpatches.Patch(facecolor=family_colors['EuSte'], edgecolor='black', label='EuSte Domains'),
        mpatches.Patch(facecolor=family_colors['HetSte'], edgecolor='black', label='HetSte Domains')
    ]
    # Shifted X from 1.02 to 1.15 to ensure it doesn't overlap the colorbar
    ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1.15, 0.5), 
              fontsize=12, frameon=False, title="Stellate Variants", title_fontsize=14)

# ==========================================
# 3. HEATMAP ENGINE
# ==========================================
def plot_triangle_heatmap(strain, df_coords):
    mat_file = f"Stellate_Identity_{strain}_vs_{strain}.tsv"
    if not os.path.exists(mat_file):
        print(f"Skipping {strain}: Missing matrix file {mat_file}")
        return
        
    print(f"Plotting {strain}...")
    df_mat = pd.read_csv(mat_file, sep='\t', index_col=0)
    matrix = df_mat.values
    N = matrix.shape[0]
    
    fig, ax = plt.subplots(figsize=(16, 10))
    fig.patch.set_facecolor('white')

    i_grid, j_grid = np.meshgrid(np.arange(N+1), np.arange(N+1), indexing='ij')
    X = (i_grid + j_grid) / 2.0
    Y = (j_grid - i_grid) / 2.0

    C = np.where(i_grid[:-1, :-1] <= j_grid[:-1, :-1], matrix, np.nan)
    mesh = ax.pcolormesh(X, Y, C, cmap='turbo', vmin=min_identity_scale, vmax=100, edgecolors='white', linewidth=0.2)

    cbar = plt.colorbar(mesh, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label('Sequence Identity (%)', rotation=270, labelpad=25, fontsize=14, fontweight='bold')
    cbar.ax.tick_params(labelsize=12)
    ax.set_aspect('equal')

    strain_coords = df_coords[df_coords['Strain'] == strain].reset_index(drop=True)
    h_array = N * 0.02
    y_array = -N * 0.025
    
    ax.plot([0, N/2], [0, N/2], color='black', linestyle='--', linewidth=1)
    ax.plot([N/2, N], [N/2, 0], color='black', linestyle='--', linewidth=1)

    for idx, row in strain_coords.iterrows():
        start = row['Matrix_Start_Index']
        end = row['Matrix_End_Index'] + 1 
        region = row['Region']
        
        btype = 'EuSte' if 'EuSte' in region else 'HetSte'
        color = family_colors[btype]
        
        # 1. Draw the actual colored block
        rect = mpatches.Rectangle((start, y_array), width=(end-start), height=h_array, 
                                  facecolor=color, edgecolor='black', lw=0.8)
        ax.add_patch(rect)
        
        mid_x = start + (end - start)/2.0
        clean_label = region.replace('EuSte ', '').replace('HetSte ', '')
        
        # Omit text entirely for visually microscopic blocks on the massive A4 plot
        if (end - start) > (N * 0.015) or strain != 'A4':
            # 2. Calculate a staggered height for the text to prevent overlapping
            # We use modulo 3 to create 3 alternating drop-down levels
            stagger_level = idx % 3
            text_y = y_array - (N * 0.02) - (stagger_level * N * 0.04)
            
            # 3. Draw the callout line extending from the block to the text
            ax.plot([mid_x, mid_x], [y_array, text_y], color='black', lw=1.0, zorder=1)
            
            # Add a small dot where the line touches the block for visual polish
            ax.plot(mid_x, y_array, marker='o', markersize=3, color='black', zorder=2)
            
            # 4. Draw horizontal text centered at the bottom of the line
            ax.text(mid_x, text_y - (N * 0.005), clean_label, ha='center', va='top', rotation=0, 
                    fontsize=10, fontweight='bold', color='#333333', 
                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=1), zorder=3)

    ax.set_title(f"{strain} X-Linked Stellate Array", fontsize=18, fontweight='bold', pad=20)
    ax.axis('off') 
    
    # Expanded the bottom limit from -0.18 to -0.25 to make room for the new callout lines
    ax.set_ylim(bottom=-N*0.25, top=N*0.65)
    ax.set_xlim(left=-N*0.05, right=N*1.05)
    
    add_figure_legend(ax)

    plt.tight_layout()
    plt.savefig(f"Final_Triangle_Stellate_{strain}.png", dpi=600, bbox_inches='tight')
    plt.savefig(f"Final_Triangle_Stellate_{strain}.svg", format='svg', dpi=600, bbox_inches='tight')
    plt.close()

# ==========================================
# 4. EXECUTION
# ==========================================
print("Generating Final Callout Heatmaps...")
for strain in strains:
    plot_triangle_heatmap(strain, df_coords)
    
print("Success!")
