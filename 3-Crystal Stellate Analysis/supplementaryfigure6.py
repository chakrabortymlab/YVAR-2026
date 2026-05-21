import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# ==========================================
# 1. LOAD STRICT DEPTH DATA
# ==========================================
male_df = pd.read_csv('male_CGcontig_depth.txt', sep='\t', header=None, names=['Contig', 'Position', 'Depth'])
female_df = pd.read_csv('female_CGcontig_depth.txt', sep='\t', header=None, names=['Contig', 'Position', 'Depth'])

CONTIG_LENGTH = 1144

# ==========================================
# 2. APPLY SMOOTHING & LOG-SCALE CLIPPING
# ==========================================
WINDOW_SIZE = 15
male_df['Smoothed'] = male_df['Depth'].rolling(window=WINDOW_SIZE, center=True, min_periods=1).mean()
female_df['Smoothed'] = female_df['Depth'].rolling(window=WINDOW_SIZE, center=True, min_periods=1).mean()

# To plot log scale properly, we must clip '0' values to a small positive baseline
MIN_LOG = 0.5
male_df['Smoothed_Log'] = male_df['Smoothed'].clip(lower=MIN_LOG)
female_df['Smoothed_Log'] = female_df['Smoothed'].clip(lower=MIN_LOG)

# ==========================================
# 3. DEFINE PRECISE COORDINATES
# ==========================================
cg40635_blocks = [(690, 1010, 98.4)]

rm_data = {
    "Copia": [(282, 405)],
    "Simple/Low Complexity": [(13, 42), (190, 221)],
    "Uncharacterized": [(0, 12), (43, 189), (222, 281), (406, 689), (1011, CONTIG_LENGTH)]
}

# Color palette
color_cg40635 = "#a50026"  
color_copia   = "#f46d43"  
color_simple  = "#abd9e9"  
color_unknown = "#e0e0e0"  
color_male    = "#1f77b4"  
color_female  = "#d62728"  

# ==========================================
# 4. SETUP FIGURE
# ==========================================
fig, (ax_cov, ax_anno) = plt.subplots(2, 1, figsize=(12, 8), gridspec_kw={'height_ratios': [2, 1]}, sharex=True)

# ==========================================
# 5. TOP TRACK: LOG-SCALE COVERAGE
# ==========================================
# Fill area under the curve using the MIN_LOG baseline
ax_cov.fill_between(male_df['Position'], MIN_LOG, male_df['Smoothed_Log'], color=color_male, alpha=0.7, label='Male Coverage (Y-linked)')
ax_cov.fill_between(female_df['Position'], MIN_LOG, female_df['Smoothed_Log'], color=color_female, alpha=0.7, label='Female Coverage (Cross-mapping)')

ax_cov.plot(male_df['Position'], male_df['Smoothed_Log'], color='#0f5282', linewidth=1.5)
ax_cov.plot(female_df['Position'], female_df['Smoothed_Log'], color='#9c1819', linewidth=1.5)

# Set Y-axis to Logarithmic
ax_cov.set_yscale('log')
ax_cov.set_ylabel("Read Coverage Depth (Log10 X)", fontsize=12, fontweight='bold')
ax_cov.set_title("Coverage and Structural Annotation of the Unmapped Contig (211000022278760)", fontsize=16, fontweight='bold', pad=15)

# Adjust Legend and Grid
ax_cov.legend(loc='upper left', fontsize=11, framealpha=0.9) # Moved to left so it doesn't block the male peak
ax_cov.grid(True, linestyle='--', alpha=0.5, which="both") # which="both" adds minor grid lines for the log scale
ax_cov.set_xlim(0, CONTIG_LENGTH)
ax_cov.set_ylim(MIN_LOG, 2000)

# ==========================================
# 6. BOTTOM TRACK: ANNOTATIONS
# ==========================================
ax_anno.set_ylim(-1.5, 3.5)
ax_anno.axis('off')

# X-axis scale
ax_anno.plot([0, CONTIG_LENGTH], [-1.0, -1.0], color='black', lw=1.2)
for tick in [0, 250, 500, 750, 1000, CONTIG_LENGTH]:
    ax_anno.plot([tick, tick], [-1.0, -1.1], color='black', lw=1.2)
    ax_anno.text(tick, -1.3, f"{tick} bp", ha='center', va='top', fontsize=10)

# FlyBase Transcript Track
ax_anno.text(-20, 2.3, "FlyBase Transcript\n(CG40635-RB)", va='center', ha='right', fontsize=11, fontweight='bold')
for (start, end, ident) in cg40635_blocks:
    rect = patches.Rectangle((start, 2.0), end - start, 0.6, 
                             linewidth=1, edgecolor='black', facecolor=color_cg40635)
    ax_anno.add_patch(rect)
    ax_anno.text(start + (end-start)/2, 2.7, f"{ident}% Seq Identity", 
                 ha='center', va='bottom', fontsize=10, color=color_cg40635, fontweight='bold')

# RepeatMasker Track
ax_anno.text(-20, 0.5, "RepeatMasker\nLandscape", va='center', ha='right', fontsize=11, fontweight='bold')

def draw_rm_blocks(regions, color, label, label_y, line_y_top):
    for i, (start, end) in enumerate(regions):
        rect = patches.Rectangle((start, 0.2), end - start, 0.6, 
                                 linewidth=1, edgecolor='black', facecolor=color)
        ax_anno.add_patch(rect)
        if i == 0 and label != "":
            midpoint = start + (end-start)/2
            ax_anno.plot([midpoint, midpoint], [0.8, line_y_top], color='black', lw=0.8)
            ax_anno.text(midpoint, label_y, label, ha='center', va='bottom', fontsize=10, fontweight='bold', color='black')

draw_rm_blocks(rm_data["Uncharacterized"], color_unknown, "", 0, 0)
draw_rm_blocks(rm_data["Simple/Low Complexity"], color_simple, "Simple/A-rich", 1.2, 0.8)
draw_rm_blocks(rm_data["Copia"], color_copia, "Copia TE", 0.9, 0.8)

# Vertical alignment guides
ax_anno.axvline(x=690, ymin=0, ymax=1, color='gray', linestyle='--', alpha=0.4)
ax_cov.axvline(x=690, ymin=0, ymax=1, color='gray', linestyle='--', alpha=0.4)

# ==========================================
# 7. FORMAT & EXPORT
# ==========================================
plt.tight_layout()
plt.savefig("CGcontig_Annotated_Coverage_LogScale.png", dpi=600, bbox_inches='tight')
plt.savefig("CGcontig_Annotated_Coverage_LogScale.svg", format='svg', bbox_inches='tight')
plt.close()

print("High-resolution, smoothed log-scale figures generated successfully.")
