import pandas as pd
import numpy as np
import os

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

df_coords = pd.DataFrame(coords, columns=["Strain", "Region", "Unit_Count", "Start", "End"])

# ==========================================
# CONFIGURATION
# ==========================================
strains = ['ISO1', 'A3', 'A4']
inter_pairs = [('A4', 'ISO1'), ('A4', 'A3'), ('ISO1', 'A3')]

# ==========================================
# 2. MATRIX CALCULATION ENGINE
# ==========================================
def get_intra_identity(mat, start, end):
    """Calculates homogenization strictly WITHIN a single array (ignoring self-hits)."""
    sub_mat = mat[start:end+1, start:end+1]
    
    # If the domain is a singleton (1 unit), naturally it is 100% identical to itself.
    if sub_mat.shape[0] < 2:
        return 100.0 
        
    # Extract the upper triangle (excluding diagonal self-100% matches)
    intra_vals = sub_mat[np.triu_indices(sub_mat.shape[0], k=1)]
    return np.mean(intra_vals)

def get_inter_identity(mat, start1, end1, start2, end2):
    """Calculates sequence identity BETWEEN two separate arrays."""
    sub_mat = mat[start1:end1+1, start2:end2+1]
    if sub_mat.size == 0:
        return np.nan
    return np.mean(sub_mat)

def get_interstrain_matrix(s1, s2):
    """Safely loads cross-strain matrix, transposing if needed."""
    f1 = f"Stellate_Identity_{s1}_vs_{s2}.tsv"
    f2 = f"Stellate_Identity_{s2}_vs_{s1}.tsv"
    if os.path.exists(f1): return pd.read_csv(f1, sep='\t', index_col=0).values
    if os.path.exists(f2): return pd.read_csv(f2, sep='\t', index_col=0).values.T
    return None

# ==========================================
# 3. EXECUTION LOGIC
# ==========================================
intra_array_stats = []
inter_array_stats = []
inter_strain_stats = []

# --- A. INTRA-STRAIN ANALYSIS ---
for strain in strains:
    mat_file = f"Stellate_Identity_{strain}_vs_{strain}.tsv"
    if not os.path.exists(mat_file): continue
    
    print(f"Processing Intra-Strain Matrix: {strain}")
    mat = pd.read_csv(mat_file, sep='\t', index_col=0).values
    strain_data = df_coords[df_coords['Strain'] == strain].to_dict('records')
    
    # 1. Homogenization WITHIN Arrays (Now captures singletons!)
    for row in strain_data:
        reg_name = f"{strain}_{row['Region'].replace(' ', '_')}"
        mean_id = get_intra_identity(mat, row['Start'], row['End'])
        
        if not np.isnan(mean_id):
            intra_array_stats.append({
                'Strain': strain,
                'Biological_Region': reg_name,
                'Unit_Count': row['Unit_Count'],
                'Mean_Internal_Identity_Pct': round(mean_id, 3)
            })
            
    # 2. Homogenization BETWEEN Arrays (Every possible intra-strain pair)
    for i in range(len(strain_data)):
        for j in range(i + 1, len(strain_data)):
            r1, r2 = strain_data[i], strain_data[j]
            reg1 = f"{strain}_{r1['Region'].replace(' ', '_')}"
            reg2 = f"{strain}_{r2['Region'].replace(' ', '_')}"
            
            mean_id = get_inter_identity(mat, r1['Start'], r1['End'], r2['Start'], r2['End'])
            if not np.isnan(mean_id):
                inter_array_stats.append({
                    'Strain': strain,
                    'Region_1': reg1,
                    'Region_2': reg2,
                    'Mean_Identity_Between_Arrays_Pct': round(mean_id, 3)
                })

# --- B. INTER-STRAIN ANALYSIS (All Combinations) ---
for s1, s2 in inter_pairs:
    mat = get_interstrain_matrix(s1, s2)
    if mat is None: continue
    
    print(f"Processing Inter-Strain Matrix: {s1} vs {s2}")
    s1_data = df_coords[df_coords['Strain'] == s1].to_dict('records')
    s2_data = df_coords[df_coords['Strain'] == s2].to_dict('records')
    
    for r1 in s1_data:
        for r2 in s2_data:
            reg1 = f"{s1}_{r1['Region'].replace(' ', '_')}"
            reg2 = f"{s2}_{r2['Region'].replace(' ', '_')}"
            
            mean_id = get_inter_identity(mat, r1['Start'], r1['End'], r2['Start'], r2['End'])
            if not np.isnan(mean_id):
                inter_strain_stats.append({
                    'Comparison': f"{s1} vs {s2}",
                    'Strain_1_Region': reg1,
                    'Strain_2_Region': reg2,
                    'Mean_InterStrain_Identity_Pct': round(mean_id, 3)
                })

# ==========================================
# 4. EXPORT RESULTS
# ==========================================
print("\nExporting Homogenization Data...")
df_intra = pd.DataFrame(intra_array_stats)
df_inter = pd.DataFrame(inter_array_stats)
df_cross = pd.DataFrame(inter_strain_stats)

df_intra.to_csv("Homogenization_Within_Arrays.tsv", sep='\t', index=False)
df_inter.to_csv("Homogenization_Between_Arrays_SameStrain.tsv", sep='\t', index=False)
df_cross.to_csv("Homogenization_AllCombinations_InterStrain.tsv", sep='\t', index=False)

print(f"\nSuccess! Exported statistics:")
print(f" - {len(df_intra)} Internal Array records (including singletons at 100%)")
print(f" - {len(df_inter)} Intra-Strain cross-combinations")
print(f" - {len(df_cross)} Inter-Strain cross-combinations")
