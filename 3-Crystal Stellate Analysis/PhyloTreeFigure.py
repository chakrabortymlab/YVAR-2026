import pandas as pd
from Bio import Phylo
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import re

# 1. Load the Data
tsv_file = "../Island_Summary_ISO-1.tsv"
tree_file = "ISO1_ML_Tree.treefile"

print("Loading data...")
islands_df = pd.read_csv(tsv_file, sep="\t")

# 2. 
color_map = {
    "SuSte_1": "#0000FF",       
    "SuSte_3": "#1E90FF",       
    "SuSte_5": "#00BFFF",       
    "SuSte_2": "#FF0000",       
    "SuSte_4": "#B22222",       
    "SuSte_6": "#8B0000",       
    "SuSte_7": "#FFA500",       
    "SuSte_8_Inverted": "#FFD700", 
    "PCKR_Forward": "#008000",   
    "PCKR_Inverted": "#32CD32"   
}

legend_names = {
    "SuSte_1": "SuSte Island 1",
    "SuSte_3": "SuSte Island 3",
    "SuSte_5": "SuSte Island 5",
    "SuSte_2": "SuSte Island 2",
    "SuSte_4": "SuSte Island 4",
    "SuSte_6": "SuSte Island 6",
    "SuSte_7": "SuSte Island 7",
    "SuSte_8_Inverted": "SuSte Island 8",
    "PCKR_Forward": "PCKR Forward",
    "PCKR_Inverted": "PCKR Inverted"
}

def get_island(unit_type, coord_start, coord_end):
    midpoint = (int(coord_start) + int(coord_end)) / 2
    for index, row in islands_df.iterrows():
        if unit_type in row['Structural_Domain']:
            if row['Start_Coordinate'] <= midpoint <= row['End_Coordinate']:
                return row['Structural_Domain']
    return "Unknown"

# 3. Parse Tree and Convert to Cladogram
print("Parsing tree...")
tree = Phylo.read(tree_file, "newick")

# Ignore mutational distances to prevent horizontal stretching
for clade in tree.find_clades():
    clade.branch_length = 1

domain_nodes = {}
for clade in tree.get_terminals():
    if clade.name:
        match = re.search(r'(PCKR|SuSte).*?_(\d+)-(\d+)$', clade.name)
        if match:
            unit_type = match.group(1)
            start = int(match.group(2))
            end = int(match.group(3))
            
            domain = get_island(unit_type, start, end)
            if domain not in domain_nodes:
                domain_nodes[domain] = []
            domain_nodes[domain].append((start, clade, domain))

print("Applying styles...")
for domain, nodes in domain_nodes.items():
    nodes.sort(key=lambda x: x[0])  
    for i, (start, clade, dom) in enumerate(nodes):
        unit_num = i + 1  
        clade.name = f"[{dom}] {unit_num}"
        clade.color = color_map.get(dom, "#000000")


fig_width = 8.5 
num_leaves = len(tree.get_terminals())
fig_height = max(11.0, num_leaves * 0.04) # Tightly spaced vertically

# Use ultra-fine lines to prevent the graphic from looking like a smudge
plt.rc('lines', linewidth=0.5)

fig, ax = plt.subplots(figsize=(fig_width, fig_height))  
Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda c: c.name if c.is_terminal() else "")

# 5. Apply Delicate Dots and Tiny Text
for text_obj in ax.texts:
    label = text_obj.get_text()
    
    # Tiny font size so 400 branches don't collide
    text_obj.set_fontsize(4) 
    
    match = re.search(r'\[(.*?)\]\s*(.*)', label)
    if match:
        domain = match.group(1)
        unit_num = match.group(2)
        
        if domain in color_map:
            color = color_map[domain]
            text_obj.set_color(color)
            
            x, y = text_obj.get_position()
            # Tiny delicate dot
            ax.plot(x, y, 'o', color=color, markersize=1.5)
            
        text_obj.set_text(f"  {unit_num}")

# Strip formatting
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_xticks([])
ax.set_yticks([])
ax.set_xlabel("")
ax.set_ylabel("")

# 6. Add Clean Legend
legend_patches = [mpatches.Patch(color=color, label=legend_names[domain]) for domain, color in color_map.items()]
plt.legend(handles=legend_patches, loc='upper left', bbox_to_anchor=(1.02, 1), 
           title="Structural Domains", fontsize=8, title_fontsize=10)

plt.subplots_adjust(left=0.01, right=0.80)

# 7. Render at 600 DPI for Print/Doc Quality
print("Saving High-Res Google Doc PNG...")
plt.savefig("ISO1_phylo.png", dpi=600, bbox_inches="tight")
plt.savefig("ISO1_phylo.svg", format="svg", bbox_inches="tight")

print("Done!")
