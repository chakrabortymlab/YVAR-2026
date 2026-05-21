import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os

def parse_rm_out(rm_file, assembled_array_size):
    """Parses a RepeatMasker .out file to create a boolean TE mask."""
    te_mask = np.zeros(assembled_array_size + 1, dtype=bool)
    if os.path.exists(rm_file):
        try:
            with open(rm_file, 'r') as f:
                for _ in range(3): next(f, None)
                for line in f:
                    parts = line.split()
                    if len(parts) >= 15:
                        start_bp, end_bp = int(parts[5]), int(parts[6])
                        s_idx = max(0, start_bp - 1)
                        e_idx = min(assembled_array_size, end_bp)
                        te_mask[s_idx:e_idx] = True
        except Exception as e:
            print(f"Warning: Error parsing {rm_file}: {e}")
    return te_mask

def create_depth_map(depth_df, assembled_array_size):
    """Creates a fast NumPy array for base-pair level math."""
    depth_map = np.zeros(assembled_array_size + 1)
    depth_map[depth_df['Pos'].values] = depth_df['Depth'].values
    return depth_map

def run_robust_unit_analysis(depth_map, bed_df, te_mask, assembled_array_size):
    results = []
    
    for idx, row in bed_df.iterrows():
        s_idx = max(0, int(row['Array_Start']) - 1)
        e_idx = min(assembled_array_size, int(row['Array_End']))
        unit_len = e_idx - s_idx
        
        unit_depths = depth_map[s_idx:e_idx]
        unit_te_mask = te_mask[s_idx:e_idx]
        valid_depths = unit_depths[(~unit_te_mask) & (unit_depths > 0)]
        
        if len(valid_depths) > 0:
            results.append({
                'Unit_ID': f"{row['Name']}_{idx}",
                'Type': 'Su(Ste)' if 'SuSte' in row['Name'] else 'PCKR',
                'Length': unit_len,
                'Median_Depth': np.median(valid_depths)
            })

    results_df = pd.DataFrame(results)
    if results_df.empty: return None, 0, 0

    full_length_suste = results_df[(results_df['Type'] == 'Su(Ste)') & (results_df['Length'] >= 2700) & (results_df['Length'] <= 2900)]
    suste_1x = full_length_suste['Median_Depth'].median() if len(full_length_suste) > 0 else 30.0
    
    full_length_pckr = results_df[(results_df['Type'] == 'PCKR') & (results_df['Length'] >= 2450) & (results_df['Length'] <= 2650)]
    pckr_1x = full_length_pckr['Median_Depth'].median() if len(full_length_pckr) > 0 else 30.0

    def calc_copies(row):
        baseline = suste_1x if row['Type'] == 'Su(Ste)' else pckr_1x
        return row['Median_Depth'] / baseline if baseline > 0 else 0

    results_df['Estimated_Copies'] = results_df.apply(calc_copies, axis=1)

    assembled_suste_bp = results_df[results_df['Type'] == 'Su(Ste)']['Length'].sum()
    assembled_pckr_bp = results_df[results_df['Type'] == 'PCKR']['Length'].sum()
    
    expected_suste_copies = results_df[results_df['Type'] == 'Su(Ste)']['Estimated_Copies'].sum()
    expected_pckr_copies = results_df[results_df['Type'] == 'PCKR']['Estimated_Copies'].sum()
    
    stats = {
        'suste_assembled_copies': len(results_df[results_df['Type'] == 'Su(Ste)']),
        'pckr_assembled_copies': len(results_df[results_df['Type'] == 'PCKR']),
        'suste_assembled_bp': assembled_suste_bp,
        'pckr_assembled_bp': assembled_pckr_bp,
        'total_assembled_bp': assembled_suste_bp + assembled_pckr_bp,
        
        'suste_expected_copies': expected_suste_copies,
        'pckr_expected_copies': expected_pckr_copies,
        'suste_expected_bp': expected_suste_copies * 2804,
        'pckr_expected_bp': expected_pckr_copies * 2584,
    }
    stats['total_expected_bp'] = stats['suste_expected_bp'] + stats['pckr_expected_bp']
    
    stats['suste_collapse'] = stats['suste_expected_copies'] / stats['suste_assembled_copies'] if stats['suste_assembled_copies'] > 0 else 0
    stats['pckr_collapse'] = stats['pckr_expected_copies'] / stats['pckr_assembled_copies'] if stats['pckr_assembled_copies'] > 0 else 0
    stats['total_collapse'] = stats['total_expected_bp'] / stats['total_assembled_bp'] if stats['total_assembled_bp'] > 0 else 0

    return stats, suste_1x, pckr_1x

def main():
    strains_config = [
        {
            "name": "ISO1",
            "hifi_depth": "../ISO1_raw_depth.tsv",  
            "illumina_depth": "ISO1_illumina_array_masked_depth.tsv",
            "bed_file": "../ISO-1.all_hits.bed",
            "gap_file": "../ISO1_gap_report.tsv",
            "rm_file": "../ISO1_Array.fasta.out", 
            "offset": 2852210
        },
        {
            "name": "A3",
            "hifi_depth": "../A3_raw_depth.tsv",
            "illumina_depth": "A3_illumina_array_masked_depth.tsv",
            "bed_file": "../A3.all_hits.bed",
            "gap_file": "../A3_gap_report.tsv",
            "rm_file": "../A3_Array.fasta.out",
            "offset": 2006621
        },
        {
            "name": "A4",
            "hifi_depth": "../A4_raw_depth.tsv",
            "illumina_depth": "A4_illumina_array_masked_depth.tsv",
            "bed_file": "../A4.all_hits.bed",
            "gap_file": "../A4_gap_report.tsv",
            "rm_file": "../A4_Array.fasta.out",
            "offset": 1975002
        }
    ]

    out_dir = "stacked_automated_figures"
    os.makedirs(out_dir, exist_ok=True)
    print(f"Outputs will be saved to: ./{out_dir}/")

    COLOR_UNIQUE = '#009E73'  
    COLOR_TE = '#E69F00'      
    COLOR_SUSTE = '#0072B2'   
    COLOR_PCKR = '#D55E00'    

    for config in strains_config:
        strain = config["name"]
        
        print(f"\n=========================================================")
        print(f"[{strain}] Starting Robust Multi-Platform Validation...")

        if not os.path.exists(config["hifi_depth"]) or not os.path.exists(config["bed_file"]):
            print(f"[{strain}] Warning: Missing required files. Skipping...")
            continue

        hifi_df = pd.read_csv(config["hifi_depth"], sep='\t', header=None, names=['Contig', 'Pos', 'Depth'])
        illumina_df = pd.read_csv(config["illumina_depth"], sep='\t', header=None, names=['Contig', 'Pos', 'Depth'])
        
        assembled_array_size = hifi_df['Pos'].max()
        hifi_map = create_depth_map(hifi_df, assembled_array_size)
        illumina_map = create_depth_map(illumina_df, assembled_array_size)
        
        bed_cols = ['Chr', 'Start', 'End', 'Name', 'Score', 'Strand', 'Identity']
        bed_df = pd.read_csv(config["bed_file"], sep='\t', header=None, names=bed_cols)
        gap_df = pd.read_csv(config["gap_file"], sep='\t') if os.path.exists(config["gap_file"]) else pd.DataFrame()
        te_mask = parse_rm_out(config["rm_file"], assembled_array_size)

        offset = config["offset"]
        array_end_coord = offset + assembled_array_size
        bed_df = bed_df[(bed_df['Start'] >= offset) & (bed_df['End'] <= array_end_coord)].copy()
        bed_df['Array_Start'] = bed_df['Start'] - offset + 1
        bed_df['Array_End'] = bed_df['End'] - offset + 1

        unit_mask = np.zeros(assembled_array_size + 1, dtype=bool)
        for _, row in bed_df.iterrows():
            s_idx = max(0, int(row['Array_Start']) - 1)
            e_idx = min(assembled_array_size, int(row['Array_End']))
            unit_mask[s_idx:e_idx] = True

        print(f"[{strain}] Running unit-by-unit collapse math...")
        hifi_stats, h_suste_1x, h_pckr_1x = run_robust_unit_analysis(hifi_map, bed_df, te_mask, assembled_array_size)
        ill_stats, i_suste_1x, i_pckr_1x = run_robust_unit_analysis(illumina_map, bed_df, te_mask, assembled_array_size)

        report_file = os.path.join(out_dir, f"{strain}_robust_collapse_report.txt")
        with open(report_file, 'w') as f:
            f.write(f"=== {strain} Y-Chromosome Megacluster Analysis ===\n\n")
            for plat_name, stats, s_1x, p_1x in [("HiFi Long-Read", hifi_stats, h_suste_1x, h_pckr_1x), 
                                                 ("Illumina Short-Read", ill_stats, i_suste_1x, i_pckr_1x)]:
                f.write(f"--- {plat_name} Platform ---\n")
                f.write(f"Baselines: Su(Ste) = {s_1x:.1f}X  |  PCKR = {p_1x:.1f}X\n\n")
                f.write(f"[ Su(Ste) Array ]\n")
                f.write(f"  Assembled: {stats['suste_assembled_copies']} copies ({stats['suste_assembled_bp']:,.0f} bp)\n")
                f.write(f"  Expected:  {stats['suste_expected_copies']:.1f} copies (~{stats['suste_expected_bp']:,.0f} bp)\n")
                f.write(f"  Collapse:  {stats['suste_collapse']:.2f}X\n\n")
                f.write(f"[ PCKR Array ]\n")
                f.write(f"  Assembled: {stats['pckr_assembled_copies']} copies ({stats['pckr_assembled_bp']:,.0f} bp)\n")
                f.write(f"  Expected:  {stats['pckr_expected_copies']:.1f} copies (~{stats['pckr_expected_bp']:,.0f} bp)\n")
                f.write(f"  Collapse:  {stats['pckr_collapse']:.2f}X\n\n")
                f.write(f"[ Total Megacluster Summary ]\n")
                f.write(f"  Total Assembled Span: {stats['total_assembled_bp']:,.0f} bp\n")
                f.write(f"  Total Expected Span:  {stats['total_expected_bp']:,.0f} bp\n")
                f.write(f"  Overall Collapse Factor: {stats['total_collapse']:.2f}X\n\n")
                
        # 4. BINNING
        bin_size_hifi = 500       
        bin_size_illumina = 2000  
        print(f"[{strain}] Binning HiFi at {bin_size_hifi}bp and Illumina at {bin_size_illumina}bp...")
        
        def bin_data(df, te_mask_array, unit_mask_array, current_bin_size):
            df_plot = df[df['Depth'] > 0].copy()
            df_plot['Window'] = (df_plot['Pos'] // current_bin_size) * current_bin_size
            binned = df_plot.groupby('Window')['Depth'].median().reset_index()
            binned['is_te'] = [np.mean(te_mask_array[w:w+current_bin_size]) > 0.50 for w in binned['Window']]
            binned['is_unit'] = [np.mean(unit_mask_array[w:w+current_bin_size]) > 0.50 for w in binned['Window']]
            return binned

        binned_hifi = bin_data(hifi_df, te_mask, unit_mask, bin_size_hifi)
        binned_illumina = bin_data(illumina_df, te_mask, unit_mask, bin_size_illumina)

        # 5. DATA PREP
        x_ill, y_ill = binned_illumina['Window'] / 1e6, binned_illumina['Depth']
        te_ill, unit_ill = binned_illumina['is_te'], binned_illumina['is_unit']
        
        x_hf, y_hf = binned_hifi['Window'] / 1e6, binned_hifi['Depth']
        te_hf, unit_hf = binned_hifi['is_te'], binned_hifi['is_unit']

        illumina_colors = np.where(~te_ill & unit_ill, COLOR_UNIQUE, COLOR_TE)
        hifi_colors = np.where(~te_hf & unit_hf, COLOR_UNIQUE, COLOR_TE)

        # ==================== ILLUMINA PLOT ====================
        fig_ill, ax_ill = plt.subplots(figsize=(24, 6))
        ax_ill.plot(x_ill, y_ill, color='dimgrey', linewidth=0.5, zorder=1, alpha=0.6)
        ax_ill.scatter(x_ill, y_ill, c=illumina_colors, s=12, zorder=2, alpha=0.9, edgecolors='none')

        ax_ill.plot([], [], marker='o', color=COLOR_UNIQUE, linestyle='None', markersize=6, label='Unique Seq')
        ax_ill.plot([], [], marker='o', color=COLOR_TE, linestyle='None', markersize=6, label='TE / Intergenic Seq')
        
        ax_ill.axhline(i_suste_1x, color=COLOR_SUSTE, linestyle='--', linewidth=2.5, alpha=0.9, label=f'Illumina Su(Ste) 1X ({i_suste_1x:.1f}X)', zorder=5)
        ax_ill.axhline(i_pckr_1x, color=COLOR_PCKR, linestyle=':', linewidth=2.5, alpha=0.9, label=f'Illumina PCKR 1X ({i_pckr_1x:.1f}X)', zorder=5)
        
        ax_ill.set_ylabel("Illumina Read Depth (X)", fontsize=14, fontweight='bold')
        ax_ill.set_xlabel("Array Position (Mbp)", fontsize=14, fontweight='bold')
        ax_ill.set_yscale('symlog', linthresh=max(min(i_suste_1x, i_pckr_1x) * 2, 10))
        max_illumina = y_ill.max()
        ax_ill.set_ylim(0, max_illumina * 1.5 if (pd.notna(max_illumina) and max_illumina > 50) else 150)

        # ====================== HIFI PLOT ======================
        fig_hf, ax_hf = plt.subplots(figsize=(24, 6))
        ax_hf.plot(x_hf, y_hf, color='dimgrey', linewidth=0.5, zorder=1, alpha=0.6)
        ax_hf.scatter(x_hf, y_hf, c=hifi_colors, s=8, zorder=2, alpha=0.9, edgecolors='none')

        ax_hf.plot([], [], marker='o', color=COLOR_UNIQUE, linestyle='None', markersize=6, label='Unique Seq')
        ax_hf.plot([], [], marker='o', color=COLOR_TE, linestyle='None', markersize=6, label='TE / Intergenic Seq')
        
        ax_hf.axhline(h_suste_1x, color=COLOR_SUSTE, linestyle='--', linewidth=2.5, alpha=0.9, label=f'HiFi Su(Ste) 1X ({h_suste_1x:.1f}X)', zorder=5)
        ax_hf.axhline(h_pckr_1x, color=COLOR_PCKR, linestyle=':', linewidth=2.5, alpha=0.9, label=f'HiFi PCKR 1X ({h_pckr_1x:.1f}X)', zorder=5)
        
        ax_hf.set_ylabel("HiFi Read Depth (X)", fontsize=14, fontweight='bold')
        ax_hf.set_xlabel("Array Position (Mbp)", fontsize=14, fontweight='bold')
        ax_hf.set_yscale('symlog', linthresh=max(min(h_suste_1x, h_pckr_1x) * 2, 10))
        max_hifi = y_hf.max()
        ax_hf.set_ylim(0, max_hifi * 1.5 if (pd.notna(max_hifi) and max_hifi > 50) else 150)

        # ================== SHARED ANNOTATIONS ==================
        for ax in (ax_ill, ax_hf):
            gap_plotted = False
            if not gap_df.empty:
                for _, row in gap_df.iterrows():
                    label = 'Assembly Gap' if not gap_plotted else ""
                    start_mb = row['Gap_Start'] / 1e6
                    end_mb = row['Gap_End'] / 1e6
                    min_vis_width = 0.01 
                    if (end_mb - start_mb) < min_vis_width:
                        midpoint = (start_mb + end_mb) / 2.0
                        start_mb = midpoint - (min_vis_width / 2.0)
                        end_mb = midpoint + (min_vis_width / 2.0)
                    ax.axvspan(start_mb, end_mb, color='yellow', alpha=0.8, lw=0, label=label, zorder=0)
                    gap_plotted = True

            pckr_plotted, suste_plotted = False, False
            for _, row in bed_df.iterrows():
                if 'SuSte' in row['Name']:
                    color, label = COLOR_SUSTE, ('Su(Ste) Unit' if not suste_plotted else "")
                    suste_plotted = True
                else:
                    color, label = COLOR_PCKR, ('PCKR Unit' if not pckr_plotted else "")
                    pckr_plotted = True
                
                unit_width = (row['Array_End'] - row['Array_Start']) / 1e6
                rect = patches.Rectangle((row['Array_Start'] / 1e6, 0), unit_width, 0.04, 
                                         color=color, transform=ax.get_xaxis_transform(), 
                                         zorder=10, label=label)
                ax.add_patch(rect)
                
            handles, labels = ax.get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            ax.legend(by_label.values(), by_label.keys(), loc='upper left', bbox_to_anchor=(1.02, 1), framealpha=0.9, fontsize=15)

        # --- SAVE PLOTS ---
        fig_ill.savefig(os.path.join(out_dir, f"{strain}_robust_illumina_coverage.png"), dpi=300, bbox_inches='tight')
        fig_ill.savefig(os.path.join(out_dir, f"{strain}_robust_illumina_coverage.svg"), format='svg', bbox_inches='tight')
        fig_hf.savefig(os.path.join(out_dir, f"{strain}_robust_hifi_coverage.png"), dpi=300, bbox_inches='tight')
        fig_hf.savefig(os.path.join(out_dir, f"{strain}_robust_hifi_coverage.svg"), format='svg', bbox_inches='tight')
        
        plt.close(fig_ill)
        plt.close(fig_hf)

if __name__ == "__main__":
    main()
