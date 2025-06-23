import streamlit as st
import matplotlib
matplotlib.use('Agg')
import pyBigWig
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import tempfile
import io
from pathlib import Path
import time
import json
from matplotlib.gridspec import GridSpec

# Set page config
st.set_page_config(
    page_title="BigWig Signal Analysis",
    page_icon="📊",
    layout="wide",
    initial_sidebar_state="expanded"
)

def save_uploaded_file(uploaded_file, temp_dir):
    """Save uploaded file to temporary directory and return path"""
    file_path = os.path.join(temp_dir, uploaded_file.name)
    with open(file_path, "wb") as f:
        f.write(uploaded_file.getbuffer())
    return file_path

def setup_custom_names(group_names, bed_names_ordered, mode="new_analysis"):
    """Setup UI for customizing BigWig and BED file names"""
    
    with st.expander("🏷️ Customize Names (Optional)", expanded=False):
        st.info("Customize display names for your files. These names will appear in plots and legends.")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("BigWig Group Names")
            custom_bigwig_names = []
            for i, name in enumerate(group_names):
                custom_name = st.text_input(
                    f"Group {i+1}:",
                    value=name,
                    key=f"{mode}_bigwig_name_{i}",
                    help=f"Original: {name}"
                )
                clean_name = custom_name.strip() if custom_name else name
                if not clean_name:
                    clean_name = name
                custom_bigwig_names.append(clean_name)
        
        with col2:
            st.subheader("BED File Names")
            custom_bed_names = []
            for i, name in enumerate(bed_names_ordered):
                custom_name = st.text_input(
                    f"BED {i+1}:",
                    value=name,
                    key=f"{mode}_bed_name_{i}",
                    help=f"Original: {name}"
                )
                clean_name = custom_name.strip() if custom_name else name
                if not clean_name:
                    clean_name = name
                custom_bed_names.append(clean_name)
        
        if len(set(custom_bigwig_names)) != len(custom_bigwig_names):
            st.warning("⚠️ BigWig names contain duplicates. This may cause plotting issues.")
        
        if len(set(custom_bed_names)) != len(custom_bed_names):
            st.warning("⚠️ BED names contain duplicates. This may cause plotting issues.")
        
        return custom_bigwig_names, custom_bed_names

def export_plot_with_format(fig, base_filename, format_type):
    """Export plot in specified format without losing the figure"""
    buf = io.BytesIO()
    
    if format_type.lower() == 'pdf':
        fig.savefig(buf, format='pdf', dpi=300, bbox_inches='tight')
        mime_type = "application/pdf"
        filename = f"{base_filename}.pdf"
    else:  # PNG
        fig.savefig(buf, format='png', dpi=300, bbox_inches='tight')
        mime_type = "image/png"
        filename = f"{base_filename}.png"
    
    buf.seek(0)
    return buf.getvalue(), filename, mime_type

def export_signal_data_to_excel(signals_data, profile_data, group_names, bed_names_ordered, 
                                bigwig_file_names, analysis_params):
    """Export all extracted signal data to Excel file with metadata"""
    
    output = io.BytesIO()
    
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        
        metadata = {
            'Parameter': [
                'Analysis Date', 'BigWig Files (Upload Order)', 'BED Files (Upload Order)',
                'BigWig Groups', 'Plot Type', 'Y-axis Maximum', 'Signal Window (bp)',
                'Max Regions per BED', 'Line Plot Extend (bp)', 'Line Plot Bin Size (bp)',
                'Heatmap Colormap', 'Heatmap Sort Regions', 'Heatmap Color Min', 'Heatmap Color Max'
            ],
            'Value': [
                pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S'),
                ' | '.join(bigwig_file_names), ' | '.join(bed_names_ordered),
                ' | '.join(group_names), analysis_params.get('plot_type', 'Unknown'),
                analysis_params.get('y_max', 'Unknown'), analysis_params.get('extend_bp', 'Unknown'),
                analysis_params.get('max_regions', 'Unknown'), analysis_params.get('line_extend', 2000),
                analysis_params.get('line_bin_size', 20),
                analysis_params.get('cmap', 'viridis'), analysis_params.get('sort_regions', True),
                analysis_params.get('vmin', 0.0), analysis_params.get('vmax', 10.0)
            ]
        }
        
        metadata_df = pd.DataFrame(metadata)
        metadata_df.to_excel(writer, sheet_name='Metadata', index=False)
        
        if signals_data:
            boxplot_data = []
            if len(group_names) == 1:
                signals_dict = signals_data[0]
                for bed_name in bed_names_ordered:
                    if bed_name in signals_dict:
                        for i, signal in enumerate(signals_dict[bed_name]):
                            boxplot_data.append({'BigWig_Group': group_names[0], 'BED_File': bed_name, 'Region_Index': i, 'Signal_Value': signal})
            else:
                for bigwig_idx, bigwig_group in enumerate(group_names):
                    signals_dict = signals_data[bigwig_idx]
                    for bed_name in bed_names_ordered:
                        if bed_name in signals_dict:
                            for i, signal in enumerate(signals_dict[bed_name]):
                                boxplot_data.append({'BigWig_Group': bigwig_group, 'BED_File': bed_name, 'Region_Index': i, 'Signal_Value': signal})
            if boxplot_data:
                pd.DataFrame(boxplot_data).to_excel(writer, sheet_name='Boxplot_Signals', index=False)

        if profile_data:
            for bigwig_idx, bigwig_group in enumerate(group_names):
                profile_dict = profile_data[bigwig_idx]
                for bed_idx, bed_name in enumerate(bed_names_ordered):
                    if bed_name in profile_dict and profile_dict[bed_name] is not None:
                        profile_info = profile_dict[bed_name]
                        base_data = {
                            'Position': profile_info['positions'], 'Mean_Signal': profile_info['mean_signal'],
                            'Std_Signal': profile_info['std_signal'], 'N_Regions': profile_info['n_regions'],
                            'BigWig_Group': bigwig_group, 'BED_File': bed_name
                        }
                        all_profiles = profile_info['all_profiles']
                        region_data = {f'Region_{i+1}': all_profiles[i, :] for i in range(min(all_profiles.shape[0], 100))}
                        all_data = {**base_data, **region_data}
                        profile_df = pd.DataFrame(all_data)
                        safe_bed_name = bed_name.replace('/', '_').replace('\\', '_').replace(' ', '_')[:20]
                        sheet_name = f'P_{bigwig_idx}_{bed_idx}_{safe_bed_name}'[:31]
                        profile_df.to_excel(writer, sheet_name=sheet_name, index=False)
    
    output.seek(0)
    return output

def load_signal_data_from_excel(uploaded_file):
    """Load signal data from uploaded Excel file"""
    
    try:
        metadata_df = pd.read_excel(uploaded_file, sheet_name='Metadata')
        metadata_dict = dict(zip(metadata_df['Parameter'], metadata_df['Value']))
        
        bigwig_file_names = metadata_dict['BigWig Files (Upload Order)'].split(' | ')
        bed_names_ordered = metadata_dict['BED Files (Upload Order)'].split(' | ')
        group_names = metadata_dict['BigWig Groups'].split(' | ')
        
        analysis_params = {
            'plot_type': metadata_dict.get('Plot Type', 'All'),
            'y_max': float(metadata_dict.get('Y-axis Maximum', 25.0)),
            'extend_bp': int(metadata_dict.get('Signal Window (bp)', 500)),
            'max_regions': int(metadata_dict.get('Max Regions per BED', 5000)),
            'line_extend': int(metadata_dict.get('Line Plot Extend (bp)', 2000)),
            'line_bin_size': int(metadata_dict.get('Line Plot Bin Size (bp)', 20)),
            'cmap': metadata_dict.get('Heatmap Colormap', 'viridis'),
            'sort_regions': bool(metadata_dict.get('Heatmap Sort Regions', True)),
            'vmin': float(metadata_dict.get('Heatmap Color Min', 0.0)),
            'vmax': float(metadata_dict.get('Heatmap Color Max', 10.0))
        }
        
        signals_data = None
        try:
            boxplot_df = pd.read_excel(uploaded_file, sheet_name='Boxplot_Signals')
            if len(group_names) == 1:
                signals_dict = {bed_name: boxplot_df[boxplot_df['BED_File'] == bed_name]['Signal_Value'].tolist() for bed_name in bed_names_ordered}
                signals_data = [signals_dict]
            else:
                signals_data = []
                for group_name in group_names:
                    group_data = boxplot_df[boxplot_df['BigWig_Group'] == group_name]
                    signals_dict = {bed_name: group_data[group_data['BED_File'] == bed_name]['Signal_Value'].tolist() for bed_name in bed_names_ordered}
                    signals_data.append(signals_dict)
        except Exception:
            signals_data = None
        
        profile_data = None
        try:
            excel_file = pd.ExcelFile(uploaded_file)
            profile_sheets = [sheet for sheet in excel_file.sheet_names if sheet.startswith('P_')]
            if profile_sheets:
                profile_data = [{} for _ in group_names]
                for sheet_name in profile_sheets:
                    try:
                        sheet_df = pd.read_excel(uploaded_file, sheet_name=sheet_name)
                        if 'BigWig_Group' in sheet_df.columns and 'BED_File' in sheet_df.columns:
                            actual_group, actual_bed = sheet_df['BigWig_Group'].iloc[0], sheet_df['BED_File'].iloc[0]
                            try:
                                group_idx = group_names.index(actual_group)
                            except ValueError: continue
                            
                            region_cols = [col for col in sheet_df.columns if col.startswith('Region_')]
                            all_profiles = np.array([sheet_df[col].values for col in region_cols]) if region_cols else np.array([sheet_df['Mean_Signal'].values])
                            
                            profile_data[group_idx][actual_bed] = {
                                'positions': sheet_df['Position'].values, 'mean_signal': sheet_df['Mean_Signal'].values,
                                'std_signal': sheet_df['Std_Signal'].values, 'all_profiles': all_profiles,
                                'n_regions': int(sheet_df['N_Regions'].iloc[0])
                            }
                    except Exception as e:
                        st.warning(f"Error reading sheet {sheet_name}: {e}")
        except Exception as e:
            st.error(f"Error loading profile data: {e}")

        return {'signals_data': signals_data, 'profile_data': profile_data, 'group_names': group_names, 
                'bed_names_ordered': bed_names_ordered, 'bigwig_file_names': bigwig_file_names, 
                'analysis_params': analysis_params, 'metadata': metadata_dict}
    
    except Exception as e:
        st.error(f"Error loading Excel file: {e}")
        return None

def setup_replicate_groups(bigwig_files):
    st.subheader("🔗 Replicate Grouping")
    st.info("Group biological replicates together. Files in the same group will be averaged.")
    if len(bigwig_files) <= 1: return [[0]] if bigwig_files else []
    if 'replicate_groups' not in st.session_state: st.session_state.replicate_groups = [[i] for i in range(len(bigwig_files))]
    st.write("**Uploaded BigWig Files:**")
    for i, bw_file in enumerate(bigwig_files): st.write(f"{i+1}. {bw_file.name}")
    col1, col2 = st.columns(2)
    if col1.button("📎 Group consecutive pairs (1+2, 3+4, etc.)"):
        st.session_state.replicate_groups = [[i, i+1] if i + 1 < len(bigwig_files) else [i] for i in range(0, len(bigwig_files), 2)]
        st.rerun()
    if col2.button("🔄 Reset to individual files"):
        st.session_state.replicate_groups = [[i] for i in range(len(bigwig_files))]
        st.rerun()
    num_groups = st.number_input("Number of groups:", 1, len(bigwig_files), len(st.session_state.replicate_groups))
    new_groups = []
    for i in range(num_groups):
        default = st.session_state.replicate_groups[i] if i < len(st.session_state.replicate_groups) else []
        selected = st.multiselect(f"Select files for Group {i+1}:", list(range(len(bigwig_files))), default, format_func=lambda x: f"{x+1}. {bigwig_files[x].name}", key=f"group_{i}")
        if selected: new_groups.append(selected)
    st.session_state.replicate_groups = new_groups
    return new_groups

def extract_signals_fast(bigwig_file_group, bed_file, extend=500, max_regions=5000):
    try:
        bed_df = pd.read_csv(bed_file, sep='\t', header=None, usecols=[0, 1, 2], names=['chr', 'start', 'end'], dtype={'chr': str, 'start': int, 'end': int})
        if len(bed_df) > max_regions: bed_df = bed_df.sample(n=max_regions, random_state=42)
    except Exception:
        return []
    try:
        bw_handles = [pyBigWig.open(p) for p in bigwig_file_group]
    except Exception:
        return []
    
    signals = []
    for _, row in bed_df.iterrows():
        center = (row['start'] + row['end']) // 2
        rep_signals = [np.mean([v for v in bw.values(row['chr'], max(0, center-extend), center+extend) if v is not None and not np.isnan(v)] or [0]) for bw in bw_handles]
        signals.append(np.mean(rep_signals) if rep_signals else 0)

    for bw in bw_handles: bw.close()
    return signals

# --- PERMANENTLY FIXED SIGNAL EXTRACTION FUNCTION ---
def extract_signals_for_profile(bigwig_file_group, bed_file, extend=2000, max_regions=5000, bin_size=20):
    """
    Extracts and bins signal for a SINGLE replicate group.
    This version uses pyBigWig.stats for robust and correct binning.
    """
    try:
        bed_df = pd.read_csv(bed_file, sep='\t', header=None, usecols=[0, 1, 2], names=['chr', 'start', 'end'], dtype={'chr': str, 'start': int, 'end': int})
        if len(bed_df) > max_regions: bed_df = bed_df.sample(n=max_regions, random_state=42)
    except Exception:
        return None

    try:
        bw_handles = [pyBigWig.open(p) for p in bigwig_file_group]
    except Exception:
        return None

    n_bins = (extend * 2) // bin_size
    all_profiles_for_all_regions = []

    for _, row in bed_df.iterrows():
        center = (row['start'] + row['end']) // 2
        start, end = max(0, center - extend), center + extend
        
        replicate_profiles = []
        for bw in bw_handles:
            try:
                # Use the more robust bw.stats for binning
                binned_vals = bw.stats(row['chr'], start, end, type="mean", nBins=n_bins, exact=False)
                # pyBigWig returns None for empty regions, replace with zeros
                replicate_profiles.append(np.nan_to_num(np.array(binned_vals, dtype=float)))
            except RuntimeError:
                # This can happen if a chromosome is not in the bigwig file
                replicate_profiles.append(np.zeros(n_bins))
        
        # Average the profiles from the replicates for this one genomic region
        if replicate_profiles:
            all_profiles_for_all_regions.append(np.mean(replicate_profiles, axis=0))

    for bw in bw_handles:
        bw.close()

    if not all_profiles_for_all_regions:
        return None

    all_profiles = np.array(all_profiles_for_all_regions)
    return {
        'positions': np.linspace(-extend, extend, n_bins),
        'mean_signal': np.mean(all_profiles, axis=0),
        'std_signal': np.std(all_profiles, axis=0),
        'all_profiles': all_profiles,
        'n_regions': len(all_profiles)
    }

def create_single_boxplot(signals_dict, bigwig_names, bed_names_ordered, y_axis_max, original_bed_names=None):
    name_mapping = dict(zip(bed_names_ordered, original_bed_names)) if original_bed_names else {n: n for n in bed_names_ordered}
    plot_data = []
    if len(bigwig_names) == 1:
        for custom_name in bed_names_ordered:
            signals = signals_dict.get(name_mapping[custom_name], [])
            plot_data.extend([{'Group': custom_name, 'Signal': s} for s in signals])
    else:
        for bw_idx, bw_name in enumerate(bigwig_names):
            for custom_name in bed_names_ordered:
                signals = signals_dict[bw_idx].get(name_mapping[custom_name], [])
                plot_data.extend([{'Group': custom_name, 'BigWig': bw_name, 'Signal': s} for s in signals])
    if not plot_data: return None
    df = pd.DataFrame(plot_data)
    fig, ax = plt.subplots(figsize=(max(8, len(bed_names_ordered) * 1.2), 6))
    if 'BigWig' in df.columns:
        sns.boxplot(data=df, x='Group', y='Signal', hue='BigWig', order=bed_names_ordered, hue_order=bigwig_names, palette='Set2', showfliers=False, ax=ax)
        ax.legend(title='BigWig Groups')
    else:
        sns.boxplot(data=df, x='Group', y='Signal', order=bed_names_ordered, hue='Group', hue_order=bed_names_ordered, palette='Set2', showfliers=False, ax=ax)
    ax.set_title('Signal Distribution at Peak Centers', fontsize=14, fontweight='bold')
    ax.set_xlabel('BED File Groups', fontsize=12); ax.set_ylabel('Mean Signal', fontsize=12)
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_ylim(0, y_axis_max)
    plt.tight_layout()
    return fig


def create_subplot_line_plot(profile_dict_list, bigwig_names, bed_names_ordered, original_bed_names=None):
    name_mapping = dict(zip(bed_names_ordered, original_bed_names)) if original_bed_names else {n: n for n in bed_names_ordered}
    valid_groups = [(c, name_mapping[c]) for c in bed_names_ordered if any(name_mapping[c] in p and p[name_mapping[c]] for p in profile_dict_list)]
    if not valid_groups: 
        return None
    
    n_groups = len(valid_groups)
    
    # Smart layout calculation
    if n_groups == 1:
        ncols, nrows = 1, 1
        fig_width, fig_height = 8, 6  # Single plot gets reasonable size
    elif n_groups == 2:
        ncols, nrows = 2, 1
        fig_width, fig_height = 12, 6  # Two plots side by side
    elif n_groups <= 4:
        ncols = min(n_groups, 2)
        nrows = (n_groups + ncols - 1) // ncols
        fig_width, fig_height = 6 * ncols, 5 * nrows
    else:
        # For many plots, use original logic but with size limits
        ncols = min(n_groups, 4)
        nrows = (n_groups + ncols - 1) // ncols
        fig_width = min(4 * ncols, 16)  # Max width of 16
        fig_height = min(4 * nrows + 1, 20)  # Max height of 20
    
    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_width, fig_height), 
                            sharey=True, squeeze=False)
    axes = axes.flatten()
    
    colors = ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D']
    
    for idx, (custom_name, original_name) in enumerate(valid_groups):
        ax = axes[idx]
        for bw_idx, bw_name in enumerate(bigwig_names):
            if original_name in profile_dict_list[bw_idx] and profile_dict_list[bw_idx][original_name]:
                p = profile_dict_list[bw_idx][original_name]
                ax.plot(p['positions'], p['mean_signal'], 
                       color=colors[bw_idx % len(colors)], 
                       label=f"{bw_name} (n={p['n_regions']})")
                if p['n_regions'] > 1:
                    sem = p['std_signal'] / np.sqrt(p['n_regions'])
                    ax.fill_between(p['positions'], p['mean_signal'] - sem, 
                                  p['mean_signal'] + sem, 
                                  color=colors[bw_idx % len(colors)], alpha=0.15)
        
        ax.set_title(custom_name, fontsize=11, fontweight='bold')
        ax.axvline(x=0, color='black', linestyle=':', alpha=0.5)
        if len(bigwig_names) > 1: 
            ax.legend(fontsize=8)
        if idx // ncols == nrows - 1: 
            ax.set_xlabel('Distance from Center (bp)')
        if idx % ncols == 0: 
            ax.set_ylabel('Mean Signal')
    
    # Hide unused subplots
    for i in range(n_groups, len(axes)): 
        axes[i].set_visible(False)
    
    # Adjust title and layout based on number of plots
    if n_groups == 1:
        fig.suptitle('Signal Profile', fontsize=14, fontweight='bold', y=0.95)
        plt.tight_layout(rect=[0, 0, 1, 0.9])
    else:
        fig.suptitle('Signal Profiles Across BED Files', fontsize=14, fontweight='bold', y=0.98)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    return fig

def create_comparison_heatmaps(profile_data, bigwig_names, bed_names, original_bed_names=None,
                               cmap='viridis', sort_regions=True, vmin=0.0, vmax=10.0):
    if not profile_data or not profile_data[0]:
        return None

    name_mapping = dict(zip(bed_names, original_bed_names)) if original_bed_names else {n: n for n in bed_names}
    
    valid_beds = [
        (c, name_mapping[c]) for c in bed_names 
        if name_mapping.get(c) in profile_data[0] and profile_data[0][name_mapping.get(c)]
    ]
    if not valid_beds:
        return None

    nrows = len(valid_beds)
    ncols = len(bigwig_names)
    
    # Smart figure sizing for heatmaps
    if ncols == 1 and nrows == 1:
        # Single heatmap
        fig_width, fig_height = 8, 6
    elif ncols == 1:
        # Single column, multiple rows
        fig_width = 8
        fig_height = min(3 * nrows, 15)  # Max height of 15
    elif nrows == 1:
        # Single row, multiple columns
        fig_width = min(5 * ncols, 20)  # Max width of 20
        fig_height = 6
    else:
        # Multiple rows and columns
        fig_width = min(4 * ncols + 1, 20)
        fig_height = min(3.5 * nrows, 18)
    
    fig = plt.figure(figsize=(fig_width, fig_height))
    
    # Adjust grid spacing based on size
    if ncols == 1:
        gs = GridSpec(nrows, ncols + 1, width_ratios=[4, 0.2], 
                     wspace=0.1, hspace=0.3)
    else:
        gs = GridSpec(nrows, ncols + 1, width_ratios=[4] * ncols + [0.2], 
                     wspace=0.15, hspace=0.4)

    im = None 

    for row_idx, (custom_bed_name, original_bed_name) in enumerate(valid_beds):
        ref_matrix = profile_data[0][original_bed_name]['all_profiles']
        if sort_regions and ref_matrix.shape[0] > 1:
            sorted_indices = np.argsort(ref_matrix.mean(axis=1))[::-1]
        else:
            sorted_indices = np.arange(ref_matrix.shape[0])

        for col_idx, bw_name in enumerate(bigwig_names):
            ax = fig.add_subplot(gs[row_idx, col_idx])

            # Adjust title font size based on number of columns
            title_fontsize = 12 if ncols <= 2 else 10
            if row_idx == 0:
                ax.set_title(bw_name, fontweight='bold', fontsize=title_fontsize)

            # Adjust ylabel font size based on number of rows
            ylabel_fontsize = 12 if nrows <= 3 else 10
            if col_idx == 0:
                ax.set_ylabel(custom_bed_name, fontweight='bold', fontsize=ylabel_fontsize)
            
            if not (col_idx < len(profile_data) and original_bed_name in profile_data[col_idx] and profile_data[col_idx][original_bed_name]):
                ax.text(0.5, 0.5, "No Data", ha='center', va='center')
                ax.set_yticks([])
                ax.set_xticks([])
            else:
                p_data = profile_data[col_idx][original_bed_name]
                matrix = p_data['all_profiles'][sorted_indices, :]
                
                im = ax.imshow(matrix, aspect='auto', interpolation='none', 
                              cmap=cmap, vmin=vmin, vmax=vmax)
                
                # Adjust annotation font size
                annotation_fontsize = 9 if ncols <= 2 else 8
                ax.text(0.02, 0.98, f"n={p_data['n_regions']}", 
                       transform=ax.transAxes, ha='left', va='top', 
                       fontsize=annotation_fontsize, 
                       bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.5))
                ax.set_yticks([])

            if row_idx == nrows - 1:
                positions = profile_data[0][original_bed_name]['positions']
                tick_pos = np.linspace(0, len(positions) - 1, 5)
                tick_labels = [f'{int(positions[int(p)]/1000)}kb' if positions[int(p)] != 0 else '0' for p in tick_pos]
                ax.set_xticks(tick_pos)
                
                # Adjust tick label font size
                tick_fontsize = 10 if ncols <= 2 else 9
                ax.set_xticklabels(tick_labels, fontsize=tick_fontsize)
                ax.set_xlabel("Distance from Center", fontsize=11)
            else:
                ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

    if im:
        cbar_ax = fig.add_subplot(gs[:, -1])
        fig.colorbar(im, cax=cbar_ax, label="Signal Intensity")

    return fig


def main():
    st.title("📊 BigWig Signal Analysis Tool")
    st.markdown("Upload BigWig/BED files or a pre-extracted Excel file to generate plots.")
    
    # Initialize session state for analysis results
    if 'analysis_results' not in st.session_state:
        st.session_state.analysis_results = None
    if 'analysis_running' not in st.session_state:
        st.session_state.analysis_running = False
    
    pre_extracted_data = None
    st.header("🔄 Import Pre-Extracted Signal Data (Optional)")
    uploaded_excel = st.file_uploader("Choose an exported signal data file (.xlsx):", type=['xlsx'], key="excel_uploader")
    
    if uploaded_excel:
        with st.spinner("Loading pre-extracted data..."):
            pre_extracted_data = load_signal_data_from_excel(uploaded_excel)
        if pre_extracted_data:
            st.success("✅ Successfully loaded pre-extracted data!")
            col1, col2 = st.columns(2)
            col1.write("**BigWig Groups:**")
            col1.json({i+1: name for i, name in enumerate(pre_extracted_data['group_names'])})
            col2.write("**BED Files:**")
            col2.json({i+1: name for i, name in enumerate(pre_extracted_data['bed_names_ordered'])})

    with st.sidebar:
        st.header("🔧 Plot Settings")
        
        plot_type_options = ["Boxplot", "Line plot", "Heatmap", "All"]
        default_index = 3 

        if pre_extracted_data:
            saved_plot_type = pre_extracted_data['analysis_params'].get('plot_type', 'All')
            if saved_plot_type == "Both": saved_plot_type = "All"
            if saved_plot_type in plot_type_options:
                default_index = plot_type_options.index(saved_plot_type)

        plot_type = st.selectbox("Select plot type:", plot_type_options, index=default_index)

        if plot_type in ["Boxplot", "All"]:
            st.subheader("Boxplot Settings")
            default_y_max = pre_extracted_data['analysis_params']['y_max'] if pre_extracted_data else 25.0
            y_max = st.number_input("Y-axis maximum:", 0.1, value=float(default_y_max))

        if plot_type in ["Heatmap", "All"]:
            st.subheader("Heatmap Settings")
            if pre_extracted_data:
                params = pre_extracted_data['analysis_params']
                default_cmap, default_sort = params.get('cmap', 'viridis'), params.get('sort_regions', True)
                default_vmin, default_vmax = params.get('vmin', 0.0), params.get('vmax', 10.0)
            else:
                default_cmap, default_sort, default_vmin, default_vmax = 'viridis', True, 0.0, 10.0
            
            cmap_options = ['viridis', 'coolwarm', 'Reds', 'Blues', 'YlGnBu', 'magma']
            cmap_index = cmap_options.index(default_cmap) if default_cmap in cmap_options else 0
            cmap_choice = st.selectbox("Colormap:", cmap_options, index=cmap_index)
            sort_regions = st.checkbox("Sort regions by mean signal (using 1st sample as reference)", value=bool(default_sort))
            vmin = st.number_input("Color scale min:", value=float(default_vmin), format="%.2f")
            vmax = st.number_input("Color scale max:", value=float(default_vmax), format="%.2f")

        st.subheader("Export Settings")
        export_format = st.selectbox("Export format:", ["PNG", "PDF"])

    if pre_extracted_data:
        st.markdown("---"); st.header("🎨 Plot Generation from Pre-Extracted Data")
        custom_bigwig_names, custom_bed_names = setup_custom_names(pre_extracted_data['group_names'], pre_extracted_data['bed_names_ordered'], mode="imported")
        
        signals_data = pre_extracted_data['signals_data']
        profile_data = pre_extracted_data['profile_data']
        original_bed_names = pre_extracted_data['bed_names_ordered']
        
        if plot_type in ["Boxplot", "All"] and signals_data:
            st.subheader("📊 Boxplot Results")
            fig = create_single_boxplot(signals_data[0] if len(custom_bigwig_names)==1 else signals_data, custom_bigwig_names, custom_bed_names, y_max, original_bed_names)
            if fig:
                st.pyplot(fig)
                f_data, f_name, f_mime = export_plot_with_format(fig, "boxplot", export_format)
                st.download_button(f"📥 Download Boxplot ({export_format})", f_data, f_name, f_mime)

        if plot_type in ["Line plot", "All"] and profile_data:
            st.subheader("📈 Line Plot Results")
            fig = create_subplot_line_plot(profile_data, custom_bigwig_names, custom_bed_names, original_bed_names)
            if fig:
                st.pyplot(fig)
                f_data, f_name, f_mime = export_plot_with_format(fig, "lineplot", export_format)
                st.download_button(f"📥 Download Line Plot ({export_format})", f_data, f_name, f_mime)
        
        if plot_type in ["Heatmap", "All"] and profile_data:
            st.subheader("🔥 Consolidated Heatmap Comparison")
            fig = create_comparison_heatmaps(profile_data, custom_bigwig_names, custom_bed_names, original_bed_names, cmap_choice, sort_regions, vmin, vmax)
            if fig:
                st.pyplot(fig)
                f_data, f_name, f_mime = export_plot_with_format(fig, "consolidated_heatmap", export_format)
                st.download_button(f"📥 Download Heatmap ({export_format})", f_data, f_name, f_mime, key=f"dl_heatmap_consolidated_imported")
        return

    st.markdown("---"); st.header("🔬 New Analysis")

    with st.sidebar:
        st.header("⚙️ Analysis Parameters")
        extend_bp = None
        if plot_type in ["Boxplot", "All"]:
             extend_bp = st.number_input("Boxplot signal window (±bp):", 50, 5000, 500, 50, disabled=st.session_state.analysis_running)
        max_regions = st.number_input("Max regions per BED:", 100, 10000, 5000, 100, disabled=st.session_state.analysis_running)

    col1, col2 = st.columns(2)
    replicate_groups = None
    with col1:
        st.subheader("📁 Upload BigWig Files")
        bigwig_files = st.file_uploader("Choose .bw files", type=['bw', 'bigwig'], accept_multiple_files=True, disabled=st.session_state.analysis_running)
        if bigwig_files and len(bigwig_files) > 1:
            if st.checkbox("🔗 Enable replicate grouping", False, disabled=st.session_state.analysis_running):
                replicate_groups = setup_replicate_groups(bigwig_files)
            else: replicate_groups = [[i] for i in range(len(bigwig_files))]
        elif bigwig_files: replicate_groups = [[0]]
    
    with col2:
        st.subheader("📄 Upload BED Files")
        bed_files = st.file_uploader("Choose .bed files", type=['bed'], accept_multiple_files=True, disabled=st.session_state.analysis_running)

    # Clear previous analysis results when new files are uploaded
    if bigwig_files or bed_files:
        if st.session_state.analysis_results:
            # Check if the uploaded files are different from the previous analysis
            current_bw_names = [f.name for f in bigwig_files] if bigwig_files else []
            current_bed_names = [f.name for f in bed_files] if bed_files else []
            
            if (current_bw_names != st.session_state.analysis_results.get('bigwig_file_names', []) or 
                current_bed_names != st.session_state.analysis_results.get('bed_file_names', [])):
                st.session_state.analysis_results = None

    if st.button("🚀 Run Analysis", type="primary", use_container_width=True, disabled=not bigwig_files or not bed_files or st.session_state.analysis_running):
        st.session_state.analysis_running = True
        
        with tempfile.TemporaryDirectory() as temp_dir:
            try:
                bigwig_paths = [save_uploaded_file(f, temp_dir) for f in bigwig_files]
                bed_paths = {Path(f.name).stem: save_uploaded_file(f, temp_dir) for f in bed_files}
                bed_names_ordered = [Path(f.name).stem for f in bed_files]
                
                grouped_bigwig_paths = [[bigwig_paths[i] for i in group] for group in replicate_groups]
                group_names = [Path(bigwig_files[g[0]].name).stem if len(g) == 1 else f"{Path(bigwig_files[g[0]].name).stem}_group" for g in replicate_groups]
                
                progress_bar = st.progress(0); status_text = st.empty()
                signals_data, profile_data = [], []
                
                total_tasks = len(grouped_bigwig_paths) * len(bed_paths)
                current_task = 0
                
                for group_idx, group_paths in enumerate(grouped_bigwig_paths):
                    signals_dict, profile_dict = {}, {}
                    for bed_name, bed_path in bed_paths.items():
                        current_task += 1
                        progress = current_task / total_tasks
                        status_text.text(f"Processing: {group_names[group_idx]} on {bed_name} ({current_task}/{total_tasks})")
                        progress_bar.progress(progress)

                        if plot_type in ["Boxplot", "All"] and extend_bp:
                            signals_dict[bed_name] = extract_signals_fast(group_paths, bed_path, extend_bp, max_regions)
                        if plot_type in ["Line plot", "Heatmap", "All"]:
                            profile_dict[bed_name] = extract_signals_for_profile(group_paths, bed_path, 2000, max_regions, 20)
                    
                    if signals_dict: signals_data.append(signals_dict)
                    if profile_dict: profile_data.append(profile_dict)
                
                status_text.text("✅ Analysis complete! Generating plots..."); progress_bar.progress(1.0)
                
                # Store results in session state
                st.session_state.analysis_results = {
                    'signals_data': signals_data,
                    'profile_data': profile_data,
                    'group_names': group_names,
                    'bed_names_ordered': bed_names_ordered,
                    'bigwig_file_names': [f.name for f in bigwig_files],
                    'bed_file_names': [f.name for f in bed_files],
                    'analysis_params': {
                        'plot_type': plot_type, 'y_max': y_max, 'extend_bp': extend_bp, 'max_regions': max_regions,
                        'line_extend': 2000, 'line_bin_size': 20, 'cmap': cmap_choice, 'sort_regions': sort_regions,
                        'vmin': vmin, 'vmax': vmax
                    }
                }

            except Exception as e:
                st.error(f"An error occurred during analysis: {e}")
                st.exception(e)
            finally:
                st.session_state.analysis_running = False

    # Display results if they exist in session state
    if st.session_state.analysis_results:
        analysis_data = st.session_state.analysis_results
        
        custom_bigwig_names, custom_bed_names = setup_custom_names(
            analysis_data['group_names'], 
            analysis_data['bed_names_ordered'], 
            mode="new_analysis"
        )
        
        if plot_type in ["Boxplot", "All"] and analysis_data['signals_data']:
            st.subheader("📊 Boxplot Results")
            fig = create_single_boxplot(
                analysis_data['signals_data'][0] if len(custom_bigwig_names)==1 else analysis_data['signals_data'], 
                custom_bigwig_names, custom_bed_names, y_max
            )
            if fig:
                st.pyplot(fig)
                f_data, f_name, f_mime = export_plot_with_format(fig, "boxplot", export_format)
                st.download_button(f"📥 Download Boxplot ({export_format})", f_data, f_name, f_mime, key="dl_boxplot_new")

        if plot_type in ["Line plot", "All"] and analysis_data['profile_data']:
            st.subheader("📈 Line Plot Results")
            fig = create_subplot_line_plot(analysis_data['profile_data'], custom_bigwig_names, custom_bed_names)
            if fig:
                st.pyplot(fig)
                f_data, f_name, f_mime = export_plot_with_format(fig, "lineplot", export_format)
                st.download_button(f"📥 Download Line Plot ({export_format})", f_data, f_name, f_mime, key="dl_lineplot_new")
        
        if plot_type in ["Heatmap", "All"] and analysis_data['profile_data']:
            st.subheader("🔥 Consolidated Heatmap Comparison")
            fig = create_comparison_heatmaps(
                analysis_data['profile_data'], custom_bigwig_names, custom_bed_names, None, 
                cmap_choice, sort_regions, vmin, vmax
            )
            if fig:
                st.pyplot(fig)
                f_data, f_name, f_mime = export_plot_with_format(fig, "consolidated_heatmap", export_format)
                st.download_button(f"📥 Download Heatmap ({export_format})", f_data, f_name, f_mime, key="dl_heatmap_consolidated_new")
        
        st.header("💾 Export Signal Data")
        excel_buffer = export_signal_data_to_excel(
            analysis_data['signals_data'], analysis_data['profile_data'], custom_bigwig_names, 
            custom_bed_names, analysis_data['bigwig_file_names'], analysis_data['analysis_params']
        )
        st.download_button(
            "📊 Download Signal Data (Excel)", excel_buffer, 
            f"signal_data_{pd.Timestamp.now().strftime('%Y%m%d_%H%M')}.xlsx", 
            "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            key="dl_excel_new"
        )

if __name__ == "__main__":
    main()