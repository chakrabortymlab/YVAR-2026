import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os

# ==========================================
# 1. PRECISE STRUCTURAL COORDINATES
# ==========================================
pckr_len = 2584
suste_len = 2806

# RepeatMasker updated coords
rm_pckr = [
    {"Repeat": "HETA", "Start": 308, "End": 421, "Color": "#4575b4"},
    {"Repeat": "Heta-1_D", "Start": 422, "End": 1631, "Color": "#74add1"},
    {"Repeat": "A-rich", "Start": 1777, "End": 1820, "Color": "#abd9e9"},
    {"Repeat": "Copia", "Start": 1867, "End": 1985, "Color": "#f46d43"}
]

# Stellate Relic removed; CG40635 Core stands alone
coding_pckr = [
    {"Region": "CG40635 Core", "Start": 2333, "End": 2584, "Color": "#a50026"}
]

rm_suste = [
    {"Repeat": "Hoppel/PROTOP", "Start": 22, "End": 1101, "Color": "#fdae61"},
    {"Repeat": "A-rich / Simple", "Start": 1472, "End": 1694, "Color": "#abd9e9"},
    {"Repeat": "Simple", "Start": 2596, "End": 2614, "Color": "#abd9e9"}
]

coding_suste = [
    {"Region": "Stellate Target", "Start": 2166, "End": 2651, "Color": "#a50026"}
]

# ==========================================
# 2. PLOTTING FUNCTION
# ==========================================
def plot_architecture(male_cov, female_cov, rm_data, coding_data, title, length, out_prefix):
    # Load data
    try:
        df_m = pd.read_csv(male_cov, sep='\t', header=None, names=['Contig', 'Position', 'Depth'])
        df_f = pd.read_csv(female_cov, sep='\t', header=None, names=['Contig', 'Position', 'Depth'])
    except Exception as e:
        print(f"Error loading coverage files: {e}")
        return

    # Baseline for log scale
    MIN_LOG = 0.5
    
    df_m['Smoothed'] = df_m['Depth'].rolling(window=15, center=True, min_periods=1).mean().clip(lower=MIN_LOG)
    df_f['Smoothed'] = df_f['Depth'].rolling(window=15, center=True, min_periods=1).mean().clip(lower=MIN_LOG)

    fig, (ax_cov, ax_anno) = plt.subplots(2, 1, figsize=(14, 8), gridspec_kw={'height_ratios': [2, 1.2]}, sharex=True)
    
    # --- TOP TRACK: COVERAGE (LOG SCALE) ---
    ax_cov.fill_between(df_m['Position'], MIN_LOG, df_m['Smoothed'], color='#1f77b4', alpha=0.7, label='Male dRNA (Y-Linked)')
    ax_cov.fill_between(df_f['Position'], MIN_LOG, df_f['Smoothed'], color='#d62728', alpha=0.7, label='Female dRNA (Cross-Mapping)')
    ax_cov.plot(df_m['Position'], df_m['Smoothed'], color='#0f5282', linewidth=1.5)
    ax_cov.plot(df_f['Position'], df_f['Smoothed'], color='#9c1819', linewidth=1.5)
    
    ax_cov.set_yscale('log')
    ax_cov.set_ylabel("Read Coverage (Log10 X)", fontsize=12, fontweight='bold')
    ax_cov.set_title(title, fontsize=16, fontweight='bold', pad=20)
    
    ax_cov.legend(loc='upper left', fontsize=11, framealpha=0.9)
    ax_cov.grid(True, linestyle='--', alpha=0.5, which="both")
    ax_cov.set_xlim(0, length)
    ax_cov.set_ylim(MIN_LOG, 2000)
    
    # --- BOTTOM TRACK: ANNOTATIONS ---
    ax_anno.set_ylim(-1.5, 3.0)
    ax_anno.axis('off')
    
    # Scale bar at y=-1.0
    ax_anno.plot([0, length], [-1.0, -1.0], color='black', lw=1.2)
    for tick in range(0, length + 1, 500):
        ax_anno.plot([tick, tick], [-1.0, -1.1], color='black', lw=1.2)
        ax_anno.text(tick, -1.3, f"{tick} bp", ha='center', va='top', fontsize=10)
        
    # Track 1: Coding Homology Domains (Top tier, y=1.5 to 2.1)
    y_code = 1.8
    h_code = 0.6
    ax_anno.text(-20, y_code, "Coding / Homology\nDomains", va='center', ha='right', fontsize=11, fontweight='bold')
    for cd in coding_data:
        rect = patches.Rectangle((cd["Start"], y_code - h_code/2), cd["End"] - cd["Start"], h_code, 
                                 linewidth=1, edgecolor='black', facecolor=cd["Color"], alpha=0.9)
        ax_anno.add_patch(rect)
        # Perfectly centered text
        ax_anno.text(cd["Start"] + (cd["End"]-cd["Start"])/2, y_code, cd["Region"], 
                     ha='center', va='center', fontsize=10, fontweight='bold', color='white')
                     
    # Track 2: RepeatMasker Landscape (Bottom tier, y=0.2 to 0.8)
    y_te = 0.5
    h_te = 0.6
    ax_anno.text(-20, y_te, "Transposable\nElements", va='center', ha='right', fontsize=11, fontweight='bold')
    for rm in rm_data:
        rect = patches.Rectangle((rm["Start"], y_te - h_te/2), rm["End"] - rm["Start"], h_te, 
                                 linewidth=1, edgecolor='black', facecolor=rm["Color"], alpha=0.9)
        ax_anno.add_patch(rect)
        
        # Stagger logic to prevent text overlap for nearby repeat elements
        v_align = 'top' if rm["Start"] % 2 == 0 else 'bottom'
        v_pos = y_te - h_te/1.7 if v_align == 'top' else y_te + h_te/1.7
        ax_anno.text(rm["Start"] + (rm["End"]-rm["Start"])/2, v_pos, rm["Repeat"], 
                     ha='center', va=v_align, fontsize=10, rotation=0)

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.15, bottom=0.15) 
    
    # Export in both formats
    plt.savefig(f"{out_prefix}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{out_prefix}.svg", format='svg', bbox_inches='tight')
    plt.close()

# ==========================================
# 3. EXECUTION
# ==========================================
print("Generating final, cleaned architectural figures...")
plot_architecture("A4_Male_PCKR_coverage.txt", "A4_Female_PCKR_coverage.txt", rm_pckr, coding_pckr,
                  "PCKR Array Transcriptional Architecture", pckr_len, "PCKR_Architecture_Final")

plot_architecture("A4_Male_SuSte_coverage.txt", "A4_Female_SuSte_coverage.txt", rm_suste, coding_suste,
                  "Su(Ste) Array Transcriptional Architecture", suste_len, "SuSte_Architecture_Final")

print("Figures Saved: PCKR_Architecture_Final.png/.svg and SuSte_Architecture_Final.png/.svg")
