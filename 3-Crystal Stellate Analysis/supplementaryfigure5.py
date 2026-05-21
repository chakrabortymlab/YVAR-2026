import matplotlib.pyplot as plt
import matplotlib.patches as patches

# ==========================================
# 1. DEFINE SHIFTED COORDINATES
# ==========================================
# The consensus is 2584 bp long. Shifted by -70 bp.
PCKR_LENGTH = 2584

# Track 1: CG40635 Transcript (Shifted to end)
cg40635_blocks = [(2263, 2583, 98.4)]

# Track 2: Release 6 Unmapped Contig (211000022278760)
# PAF mapping shifted by -70 bp (modulo 2584)
contig_blocks = [
    (1573, 2510), # Hit 1 shifted
    (2517, PCKR_LENGTH), # Hit 2 shifted (tail)
    (0, 259)             # Hit 2 shifted (wrapped head)
]

# Track 4: RepeatMasker Data (Shifted by -70 bp)
rm_data = {
    "HeT-A": [(238, 1561)], 
    "Copia": [(1797, 1915)],
    "Simple/Low Complexity": [(1707, 1750)],
    "Uncharacterized": [(0, 237), (1562, 1706), (1751, 1796), (1916, 2262)]
}

# Colors
color_cg40635 = "#a50026" # Dark red for Transcript
color_contig  = "#9970ab" # Purple for Assembly Contig
color_heta    = "#4575b4" # Blue for HeT-A
color_copia   = "#f46d43" # Orange for Copia
color_simple  = "#abd9e9" # Light blue for simple repeats
color_unknown = "#e0e0e0" # Grey for uncharacterized

# ==========================================
# 2. SETUP FIGURE
# ==========================================
fig, ax = plt.subplots(figsize=(12, 7)) # Increased height for 4 tracks
ax.set_xlim(-150, PCKR_LENGTH + 150)
ax.set_ylim(-2.5, 4.5)
ax.axis('off')

# Draw X-axis scale
ax.plot([0, PCKR_LENGTH], [-1.5, -1.5], color='black', lw=1.2)
for tick in [0, 500, 1000, 1500, 2000, PCKR_LENGTH]:
    ax.plot([tick, tick], [-1.5, -1.6], color='black', lw=1.2)
    ax.text(tick, -1.75, f"{tick} bp", ha='center', va='top', fontsize=10)

# ==========================================
# 3. TRACK 1: CG40635 Homology (Top)
# ==========================================
ax.text(-50, 2.9, "FlyBase Transcript\n(CG40635-RB)", va='center', ha='right', fontsize=11, fontweight='bold')

for (start, end, ident) in cg40635_blocks:
    rect = patches.Rectangle((start, 2.6), end - start, 0.6, 
                             linewidth=1, edgecolor='black', facecolor=color_cg40635)
    ax.add_patch(rect)
    ax.text(start + (end-start)/2, 3.4, f"{ident}% Sequence Identity", 
            ha='center', va='bottom', fontsize=10, color=color_cg40635, fontweight='bold')

# ==========================================
# 4. TRACK 2: Unmapped Release 6 Contig
# ==========================================
ax.text(-50, 1.5, "Release 6 Contig\n(211...78760)", va='center', ha='right', fontsize=11, fontweight='bold')

for (start, end) in contig_blocks:
    rect = patches.Rectangle((start, 1.2), end - start, 0.6, 
                             linewidth=1, edgecolor='black', facecolor=color_contig)
    ax.add_patch(rect)

# Draw dashed line to show contig wrapping over the arbitrary consensus boundary
ax.plot([contig_blocks[1][1], PCKR_LENGTH + 50, -50, contig_blocks[2][0]], 
        [1.5, 1.5, 1.5, 1.5], color='black', linestyle='--', lw=1.2, zorder=0)

ax.text(1573 + (2510-1573)/2, 2.0, "Truncated Assembly Fragment (~1.1 kb)", 
        ha='center', va='bottom', fontsize=10, color=color_contig, fontweight='bold')

# Vertical connecting lines from Transcript down to Contig
ax.plot([cg40635_blocks[0][0], cg40635_blocks[0][0]], [1.8, 2.6], color='black', linestyle=':', lw=1, alpha=0.5)
ax.plot([cg40635_blocks[0][1], cg40635_blocks[0][1]], [1.8, 2.6], color='black', linestyle=':', lw=1, alpha=0.5)

# ==========================================
# 5. TRACK 3: PCKR Consensus Backbone (Middle)
# ==========================================
ax.text(-50, 0, "Resolved PCKR\nConsensus Unit", va='center', ha='right', fontsize=11, fontweight='bold')
backbone = patches.Rectangle((0, -0.1), PCKR_LENGTH, 0.2, color='black', zorder=1)
ax.add_patch(backbone)

# Vertical connecting lines from Contig down to Backbone
ax.plot([contig_blocks[0][0], contig_blocks[0][0]], [0.1, 1.2], color='black', linestyle=':', lw=1, alpha=0.5)

# ==========================================
# 6. TRACK 4: RepeatMasker Annotation (Bottom)
# ==========================================
ax.text(-50, -0.8, "RepeatMasker\nLandscape", va='center', ha='right', fontsize=11, fontweight='bold')

def draw_rm_blocks(regions, color, label, label_y, line_y_top):
    for i, (start, end) in enumerate(regions):
        rect = patches.Rectangle((start, -1.1), end - start, 0.6, 
                                 linewidth=1, edgecolor='black', facecolor=color)
        ax.add_patch(rect)
        if i == 0 and label != "":
            midpoint = start + (end-start)/2
            ax.plot([midpoint, midpoint], [line_y_top, label_y + 0.1], color='black', lw=0.8)
            ax.text(midpoint, label_y, label, ha='center', va='top', fontsize=10, fontweight='bold', color='black')

draw_rm_blocks(rm_data["Uncharacterized"], color_unknown, "", 0, 0)
draw_rm_blocks(rm_data["HeT-A"], color_heta, "HeT-A ~1.3 kb", -1.2, -1.1)
draw_rm_blocks(rm_data["Simple/Low Complexity"], color_simple, "A-rich", -2.2, -1.1)
draw_rm_blocks(rm_data["Copia"], color_copia, "Copia", -1.2, -1.1)

# ==========================================
# 7. FORMAT & EXPORT
# ==========================================
plt.title("Structural Architecture of the PCKR Array Unit", fontsize=16, fontweight='bold', pad=30)
plt.tight_layout()

plt.savefig("PCKR_Architecture_With_Contig.png", dpi=600, bbox_inches='tight')
plt.savefig("PCKR_Architecture_With_Contig.svg", format='svg', bbox_inches='tight')
plt.close()

print("Figure generation complete. Includes FlyBase contig mapping.")
