import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os

# ==========================================
# CONFIGURATION
# ==========================================
strains = ['ISO-1', 'A3', 'A4']
inter_pairs = [('A4', 'ISO-1'), ('A4', 'A3'), ('ISO-1', 'A3')]

# Updated Island sizes to reflect the new 8-island configuration
# Island 6 is split to capture the ENTIRE distal transition zone (the full "dark triangle") into Island 7
suste_island_sizes = {
    'ISO-1': [40, 67, 46, 27, 52, 60, 14, 52],
    'A3': [36, 39, 48, 72, 64, 48, 15, 48],
    'A4': [39, 98, 46, 73, 53, 70, 14, 50]
}

# ==========================================
# ANNOTATION LOGIC
# ==========================================
def get_label_info(labels):
    info = []
    # Keep track of unit counts per strain to dynamically determine island boundaries
    counts = {s: 0 for s in strains}
    
    cumulative_thresholds = {}
    for s in suste_island_sizes:
        sizes = suste_island_sizes[s]
        cumulative_thresholds[s] = [sum(sizes[:i+1]) for i in range(len(sizes))]
        
    for lbl in labels:
        strain = lbl.split('_')[0]
        if 'PCKR' in lbl.upper():
            info.append({'strain': strain, 'type': 'PCKR'})
        else:
            counts[strain] += 1
            island_idx = 1
            for i, thresh in enumerate(cumulative_thresholds[strain]):
                if counts[strain] <= thresh:
                    island_idx = i + 1
                    break
            info.append({'strain': strain, 'type': island_idx})
    return info

# ==========================================
# STATS AND PLOTTING
# ==========================================
def run_all_stats_and_plot():
    all_data = []

    # 1. PROCESS INTRA-STRAIN
    for strain in strains:
        f_name = f'Matrix_Intra_{strain}.tsv'
        if not os.path.exists(f_name): continue
        
        print(f"Processing Intra-Strain: {strain}...")
        df = pd.read_csv(f_name, sep='\t', index_col=0)
        mat = df.values
        labels = list(df.index)
        info = get_label_info(labels)
        
        intra_vals, adj_vals = [], []
        N = len(labels)
        for i in range(N):
            for j in range(i+1, N):
                t1, t2 = info[i]['type'], info[j]['type']
                if t1 == 'PCKR' or t2 == 'PCKR': continue
                    
                val = mat[i, j]
                if t1 == t2:
                    intra_vals.append(val)
                    all_data.append({'Dataset': strain, 'Category': 'Aligned Domains (Intra/Orthologous)', 'Identity': val})
                elif abs(t1 - t2) == 1:
                    adj_vals.append(val)
                    all_data.append({'Dataset': strain, 'Category': 'Adjacent Boundaries (Shifted)', 'Identity': val})
                    
        u_stat, p_val = stats.mannwhitneyu(intra_vals, adj_vals, alternative='greater')
        print(f"  {strain} Mann-Whitney p-value: {p_val:.2e}\n")

    # 2. PROCESS INTER-STRAIN
    for s1, s2 in inter_pairs:
        f_name = f'Matrix_Inter_{s1}_vs_{s2}.tsv'
        if not os.path.exists(f_name): continue
        
        dataset_label = f"{s1} vs {s2}"
        print(f"Processing Inter-Strain: {dataset_label}...")
        df = pd.read_csv(f_name, sep='\t', index_col=0)
        mat = df.values
        labels = list(df.index)
        info = get_label_info(labels)
        
        ortho_vals, adj_vals = [], []
        N = len(labels)
        for i in range(N):
            for j in range(i+1, N):
                inf1, inf2 = info[i], info[j]
                
                # Only compare across the two strains
                if inf1['strain'] == inf2['strain']: continue
                    
                t1, t2 = inf1['type'], inf2['type']
                if t1 == 'PCKR' or t2 == 'PCKR': continue
                    
                val = mat[i, j]
                if t1 == t2:
                    ortho_vals.append(val)
                    all_data.append({'Dataset': dataset_label, 'Category': 'Aligned Domains (Intra/Orthologous)', 'Identity': val})
                elif abs(t1 - t2) == 1:
                    adj_vals.append(val)
                    all_data.append({'Dataset': dataset_label, 'Category': 'Adjacent Boundaries (Shifted)', 'Identity': val})
                    
        u_stat, p_val = stats.mannwhitneyu(ortho_vals, adj_vals, alternative='greater')
        print(f"  {dataset_label} Mann-Whitney p-value: {p_val:.2e}\n")

    # 3. GENERATE MASTER VIOLIN PLOT
    plot_df = pd.DataFrame(all_data)
    
    # Set the order so Intra goes first, then Inter
    order = strains + [f"{s1} vs {s2}" for s1, s2 in inter_pairs]
    
    fig, ax = plt.subplots(figsize=(14, 7))
    fig.patch.set_facecolor('white')

    sns.violinplot(
        data=plot_df, 
        x='Dataset', 
        y='Identity', 
        hue='Category', 
        split=True, 
        inner='quartile',
        order=order,
        palette={'Aligned Domains (Intra/Orthologous)': '#6baed6', 'Adjacent Boundaries (Shifted)': '#fd8d3c'},
        linewidth=1.2,
        ax=ax
    )

    # Vertical divider separating Intra from Inter
    ax.axvline(2.5, color='black', linestyle='--', alpha=0.5, linewidth=2)
    ax.text(1, 100.5, "Within-Strain (Intra-Array)", fontsize=14, fontweight='bold', ha='center', color='#333333')
    ax.text(4, 100.5, "Between-Strain (Orthologous Arrays)", fontsize=14, fontweight='bold', ha='center', color='#333333')

    ax.set_title("Statistical Distribution of Su(Ste) Sequence Homology", fontsize=18, fontweight='bold', pad=25)
    ax.set_ylabel("Pairwise Sequence Identity (%)", fontsize=14, fontweight='bold')
    ax.set_xlabel("Genomic Comparison", fontsize=14, fontweight='bold')
    
    # Zoom in on the relevant 85-100% homogenization range
    ax.set_ylim(85, 101.5) 
    
    ax.grid(axis='y', linestyle=':', alpha=0.7)
    ax.legend(title="Spatial / Structural Relationship", loc='lower center', bbox_to_anchor=(0.5, -0.2), ncol=2, frameon=False, fontsize=12, title_fontsize=12)

    plt.tight_layout()
    plt.savefig("Statistical_Distributions_Master_Violin.png", dpi=600, bbox_inches='tight')
    plt.savefig("Statistical_Distributions_Master_Violin.svg", format='svg', dpi=600, bbox_inches='tight')
    plt.close()
    print("Success! Master Violin Plot generated.")

if __name__ == "__main__":
    run_all_stats_and_plot()
