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
                # Validate custom name
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
                # Validate custom name
                clean_name = custom_name.strip() if custom_name else name
                if not clean_name:
                    clean_name = name
                custom_bed_names.append(clean_name)
        
        # Validate for duplicates
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
        
        # Sheet 1: Metadata
        metadata = {
            'Parameter': [
                'Analysis Date',
                'BigWig Files (Upload Order)', 
                'BED Files (Upload Order)',
                'BigWig Groups',
                'Plot Type',
                'Y-axis Maximum',
                'Signal Window (bp)',
                'Max Regions per BED',
                'Line Plot Extend (bp)',
                'Line Plot Bin Size (bp)'
            ],
            'Value': [
                pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S'),
                ' | '.join(bigwig_file_names),
                ' | '.join(bed_names_ordered),
                ' | '.join(group_names),
                analysis_params.get('plot_type', 'Unknown'),
                analysis_params.get('y_max', 'Unknown'),
                analysis_params.get('extend_bp', 'Unknown'),
                analysis_params.get('max_regions', 'Unknown'),
                analysis_params.get('line_extend', 2000),
                analysis_params.get('line_bin_size', 20)
            ]
        }
        
        metadata_df = pd.DataFrame(metadata)
        metadata_df.to_excel(writer, sheet_name='Metadata', index=False)
        
        # Sheet 2: Boxplot Signal Data
        if signals_data:
            boxplot_data = []
            
            if len(group_names) == 1:
                # Single BigWig group
                signals_dict = signals_data[0]
                for bed_name in bed_names_ordered:
                    if bed_name in signals_dict:
                        for i, signal in enumerate(signals_dict[bed_name]):
                            boxplot_data.append({
                                'BigWig_Group': group_names[0],
                                'BED_File': bed_name,
                                'Region_Index': i,
                                'Signal_Value': signal
                            })
            else:
                # Multiple BigWig groups
                for bigwig_idx, bigwig_group in enumerate(group_names):
                    signals_dict = signals_data[bigwig_idx]
                    for bed_name in bed_names_ordered:
                        if bed_name in signals_dict:
                            for i, signal in enumerate(signals_dict[bed_name]):
                                boxplot_data.append({
                                    'BigWig_Group': bigwig_group,
                                    'BED_File': bed_name,
                                    'Region_Index': i,
                                    'Signal_Value': signal
                                })
            
            if boxplot_data:
                boxplot_df = pd.DataFrame(boxplot_data)
                boxplot_df.to_excel(writer, sheet_name='Boxplot_Signals', index=False)
        
        # Sheet 3: Line Plot Profile Data - FIXED VERSION
        if profile_data:
            for bigwig_idx, bigwig_group in enumerate(group_names):
                profile_dict = profile_data[bigwig_idx]
                
                for bed_name in bed_names_ordered:
                    if bed_name in profile_dict and profile_dict[bed_name] is not None:
                        profile_info = profile_dict[bed_name]
                        
                        # Create base DataFrame
                        base_data = {
                            'Position': profile_info['positions'],
                            'Mean_Signal': profile_info['mean_signal'],
                            'Std_Signal': profile_info['std_signal'],
                            'N_Regions': profile_info['n_regions'],
                            'BigWig_Group': bigwig_group,
                            'BED_File': bed_name
                        }
                        
                        # FIXED: Create region columns efficiently using pd.concat
                        all_profiles = profile_info['all_profiles']
                        region_data = {}
                        for i in range(min(all_profiles.shape[0], 100)):  # Limit to first 100 regions
                            region_data[f'Region_{i+1}'] = all_profiles[i, :]
                        
                        # Combine all data at once
                        all_data = {**base_data, **region_data}
                        profile_df = pd.DataFrame(all_data)
                        
                        # Create safe sheet name
                        safe_bigwig = bigwig_group.replace('/', '_').replace('\\', '_')[:10]
                        safe_bed = bed_name.replace('/', '_').replace('\\', '_')[:15]
                        sheet_name = f'Profile_{safe_bigwig}_{safe_bed}'[:31]  # Excel sheet name limit
                        
                        profile_df.to_excel(writer, sheet_name=sheet_name, index=False)
    
    output.seek(0)
    return output

def load_signal_data_from_excel(uploaded_file):
    """Load signal data from uploaded Excel file"""
    
    try:
        # Read metadata
        metadata_df = pd.read_excel(uploaded_file, sheet_name='Metadata')
        metadata_dict = dict(zip(metadata_df['Parameter'], metadata_df['Value']))
        
        # Parse metadata
        bigwig_file_names = metadata_dict['BigWig Files (Upload Order)'].split(' | ')
        bed_names_ordered = metadata_dict['BED Files (Upload Order)'].split(' | ')
        group_names = metadata_dict['BigWig Groups'].split(' | ')
        
        analysis_params = {
            'plot_type': metadata_dict.get('Plot Type', 'Both'),
            'y_max': float(metadata_dict.get('Y-axis Maximum', 25.0)),
            'extend_bp': int(metadata_dict.get('Signal Window (bp)', 500)),
            'max_regions': int(metadata_dict.get('Max Regions per BED', 5000)),
            'line_extend': int(metadata_dict.get('Line Plot Extend (bp)', 2000)),
            'line_bin_size': int(metadata_dict.get('Line Plot Bin Size (bp)', 20))
        }
        
        # Read boxplot data
        signals_data = None
        try:
            boxplot_df = pd.read_excel(uploaded_file, sheet_name='Boxplot_Signals')
            
            if len(group_names) == 1:
                # Single BigWig group
                signals_dict = {}
                for bed_name in bed_names_ordered:
                    bed_data = boxplot_df[boxplot_df['BED_File'] == bed_name]
                    if not bed_data.empty:
                        signals_dict[bed_name] = bed_data['Signal_Value'].tolist()
                signals_data = [signals_dict]
            else:
                # Multiple BigWig groups
                signals_data = []
                for group_name in group_names:
                    signals_dict = {}
                    group_data = boxplot_df[boxplot_df['BigWig_Group'] == group_name]
                    for bed_name in bed_names_ordered:
                        bed_data = group_data[group_data['BED_File'] == bed_name]
                        if not bed_data.empty:
                            signals_dict[bed_name] = bed_data['Signal_Value'].tolist()
                    signals_data.append(signals_dict)
        
        except Exception as e:
            st.warning(f"Could not load boxplot data: {e}")
            signals_data = None
        
        # FIXED: Read profile data with better debugging
        profile_data = None
        try:
            # Get all sheet names that start with 'Profile_'
            excel_file = pd.ExcelFile(uploaded_file)
            profile_sheets = [sheet for sheet in excel_file.sheet_names if sheet.startswith('Profile_')]
            
            st.write(f"DEBUG: Found {len(profile_sheets)} profile sheets: {profile_sheets}")
            
            if profile_sheets:
                profile_data = []
                
                for group_idx, group_name in enumerate(group_names):
                    profile_dict = {}
                    
                    st.write(f"DEBUG: Processing group {group_idx}: {group_name}")
                    
                    for bed_name in bed_names_ordered:
                        st.write(f"DEBUG: Looking for BED {bed_name} in group {group_name}")
                        
                        # Try multiple matching strategies
                        matching_sheet = None
                        
                        # Strategy 1: Exact match
                        safe_bigwig = group_name.replace('/', '_').replace('\\', '_')[:10]
                        safe_bed = bed_name.replace('/', '_').replace('\\', '_')[:15]
                        exact_pattern = f'Profile_{safe_bigwig}_{safe_bed}'
                        
                        for sheet in profile_sheets:
                            if sheet == exact_pattern or sheet.startswith(exact_pattern):
                                matching_sheet = sheet
                                break
                        
                        # Strategy 2: Partial match on group name
                        if not matching_sheet:
                            for sheet in profile_sheets:
                                if safe_bigwig in sheet and safe_bed in sheet:
                                    matching_sheet = sheet
                                    break
                        
                        # Strategy 3: Look for any sheet that contains both group and bed identifiers
                        if not matching_sheet:
                            for sheet in profile_sheets:
                                # Extract parts from the original names for better matching
                                group_parts = group_name.replace('_', ' ').split()
                                bed_parts = bed_name.replace('_', ' ').split()
                                
                                # Look for meaningful parts (longer than 2 chars)
                                group_keywords = [part for part in group_parts if len(part) > 2]
                                bed_keywords = [part for part in bed_parts if len(part) > 2]
                                
                                group_match = any(keyword.lower() in sheet.lower() for keyword in group_keywords)
                                bed_match = any(keyword.lower() in sheet.lower() for keyword in bed_keywords)
                                
                                if group_match and bed_match:
                                    matching_sheet = sheet
                                    break
                        
                        if matching_sheet:
                            st.write(f"DEBUG: Found matching sheet: {matching_sheet}")
                            try:
                                profile_df = pd.read_excel(uploaded_file, sheet_name=matching_sheet)
                                
                                # Extract individual region profiles
                                region_cols = [col for col in profile_df.columns if col.startswith('Region_')]
                                all_profiles = []
                                for col in region_cols:
                                    all_profiles.append(profile_df[col].values)
                                
                                if all_profiles:
                                    all_profiles = np.array(all_profiles)
                                else:
                                    # Create dummy profiles from mean
                                    all_profiles = np.array([profile_df['Mean_Signal'].values])
                                
                                profile_info = {
                                    'positions': profile_df['Position'].values,
                                    'mean_signal': profile_df['Mean_Signal'].values,
                                    'std_signal': profile_df['Std_Signal'].values,
                                    'all_profiles': all_profiles,
                                    'n_regions': int(profile_df['N_Regions'].iloc[0])
                                }
                                
                                profile_dict[bed_name] = profile_info
                                st.write(f"DEBUG: Successfully loaded profile for {bed_name} with n_regions={profile_info['n_regions']}")
                                
                            except Exception as e:
                                st.error(f"Error loading profile data from sheet {matching_sheet}: {e}")
                        else:
                            st.warning(f"DEBUG: No matching sheet found for group {group_name}, bed {bed_name}")
                    
                    profile_data.append(profile_dict)
                    st.write(f"DEBUG: Group {group_idx} profile_dict keys: {list(profile_dict.keys())}")
        
        except Exception as e:
            st.error(f"Error loading profile data: {e}")
            profile_data = None
        
        return {
            'signals_data': signals_data,
            'profile_data': profile_data,
            'group_names': group_names,
            'bed_names_ordered': bed_names_ordered,
            'bigwig_file_names': bigwig_file_names,
            'analysis_params': analysis_params,
            'metadata': metadata_dict
        }
    
    except Exception as e:
        st.error(f"Error loading Excel file: {e}")
        return None

def setup_replicate_groups(bigwig_files):
    """Setup UI for defining replicate groups"""
    st.subheader("🔗 Replicate Grouping")
    st.info("Group biological replicates together. Files in the same group will be averaged during signal extraction.")
    
    if len(bigwig_files) <= 1:
        return [[0]] if bigwig_files else []
    
    # Initialize session state for groups
    if 'replicate_groups' not in st.session_state:
        st.session_state.replicate_groups = [[i] for i in range(len(bigwig_files))]
    
    # Display current files
    st.write("**Uploaded BigWig Files:**")
    for i, bw_file in enumerate(bigwig_files):
        st.write(f"{i+1}. {bw_file.name}")
    
    # Quick grouping options
    col1, col2 = st.columns(2)
    
    with col1:
        if st.button("📎 Group consecutive pairs (1+2, 3+4, etc.)"):
            groups = []
            for i in range(0, len(bigwig_files), 2):
                if i + 1 < len(bigwig_files):
                    groups.append([i, i + 1])
                else:
                    groups.append([i])
            st.session_state.replicate_groups = groups
            st.rerun()
    
    with col2:
        if st.button("🔄 Reset to individual files"):
            st.session_state.replicate_groups = [[i] for i in range(len(bigwig_files))]
            st.rerun()
    
    # Manual group assignment
    num_groups = st.number_input("Number of groups:", min_value=1, max_value=len(bigwig_files), 
                                value=len(st.session_state.replicate_groups))
    
    new_groups = []
    for group_idx in range(num_groups):
        st.write(f"**Group {group_idx + 1}:**")
        
        # Get current group or initialize empty
        current_group = st.session_state.replicate_groups[group_idx] if group_idx < len(st.session_state.replicate_groups) else []
        
        # Multi-select for files in this group
        selected_files = st.multiselect(
            f"Select files for Group {group_idx + 1}:",
            options=list(range(len(bigwig_files))),
            default=current_group,
            format_func=lambda x: f"{x+1}. {bigwig_files[x].name}",
            key=f"group_{group_idx}"
        )
        
        if selected_files:
            new_groups.append(selected_files)
            file_names = [bigwig_files[i].name for i in selected_files]
            if len(selected_files) > 1:
                st.success(f"✅ Will average: {', '.join(file_names)}")
            else:
                st.info(f"ℹ️ Single file: {file_names[0]}")
    
    # Validate groups
    all_assigned = set()
    for group in new_groups:
        all_assigned.update(group)
    
    unassigned = set(range(len(bigwig_files))) - all_assigned
    if unassigned:
        st.warning(f"⚠️ Unassigned files: {[bigwig_files[i].name for i in unassigned]}")
    
    # Check for duplicates
    duplicates = []
    for i, group in enumerate(new_groups):
        if len(set(group)) != len(group):
            duplicates.append(i + 1)
    
    if duplicates:
        st.error(f"❌ Groups {duplicates} contain duplicate files!")
        return None
    
    # Update session state
    st.session_state.replicate_groups = new_groups
    
    return new_groups

# [Keep all the extraction functions exactly the same as before]
def extract_signals_fast(bigwig_file_groups, bed_file, extend=500, max_regions=5000):
    """Extract signals for boxplots - peak center only, with replicate averaging"""
    
    # Read BED file
    try:
        bed_df = pd.read_csv(bed_file, sep='\t', header=None, 
                            usecols=[0, 1, 2],
                            names=['chr', 'start', 'end'],
                            dtype={'chr': str, 'start': int, 'end': int})
        
        original_count = len(bed_df)
        
        if len(bed_df) > max_regions:
            bed_df = bed_df.sample(n=max_regions, random_state=42)
        
    except Exception as e:
        st.error(f"ERROR reading BED file {bed_file}: {e}")
        return []
    
    # Open all bigwig files
    try:
        all_bw_handles = []
        for group in bigwig_file_groups:
            group_handles = [pyBigWig.open(group[i]) for i in range(len(group))]
            all_bw_handles.append(group_handles)
    except Exception as e:
        st.error(f"ERROR opening BigWig files: {e}")
        return []
    
    signals = []
    valid_regions = 0
    
    for idx, row in bed_df.iterrows():
        try:
            center = (row['start'] + row['end']) // 2
            region_start = max(0, center - extend)
            region_end = center + extend
            
            # Extract signals from each replicate group
            group_signals = []
            for group_handles in all_bw_handles:
                # Average within each replicate group
                replicate_signals = []
                for bw in group_handles:
                    values = bw.values(row['chr'], region_start, region_end)
                    if values is not None:
                        clean_values = [v for v in values if v is not None and not np.isnan(v)]
                        if clean_values:
                            mean_signal = np.mean(clean_values)
                            replicate_signals.append(mean_signal)
                        else:
                            replicate_signals.append(0)
                    else:
                        replicate_signals.append(0)
                
                # Average across replicates in this group
                if replicate_signals:
                    group_consensus = np.mean(replicate_signals)
                    group_signals.append(group_consensus)
                else:
                    group_signals.append(0)
            
            # Average across all groups (if multiple conditions)
            if group_signals:
                final_signal = np.mean(group_signals)
                signals.append(final_signal)
                valid_regions += 1
            else:
                signals.append(0)
                
        except Exception as e:
            signals.append(0)
    
    # Close all bigwig files
    for group_handles in all_bw_handles:
        for bw in group_handles:
            bw.close()
        
    return signals

def extract_signals_for_profile(bigwig_file_groups, bed_file, extend=2000, max_regions=5000, bin_size=20):
    """Extract signals for line plots - returns position-wise signals with replicate averaging"""
    
    # Read BED file
    try:
        bed_df = pd.read_csv(bed_file, sep='\t', header=None, 
                            usecols=[0, 1, 2],
                            names=['chr', 'start', 'end'],
                            dtype={'chr': str, 'start': int, 'end': int})
        
        original_count = len(bed_df)
        
        if len(bed_df) > max_regions:
            bed_df = bed_df.sample(n=max_regions, random_state=42)
        
    except Exception as e:
        st.error(f"ERROR reading BED file {bed_file}: {e}")
        return None
    
    # Open all bigwig files
    try:
        all_bw_handles = []
        for group in bigwig_file_groups:
            group_handles = [pyBigWig.open(group[i]) for i in range(len(group))]
            all_bw_handles.append(group_handles)
    except Exception as e:
        st.error(f"ERROR opening BigWig files: {e}")
        return None
    
    # Calculate number of bins
    total_length = extend * 2
    n_bins = total_length // bin_size
    
    # Store all profiles
    all_profiles = []
    valid_regions = 0
    
    for idx, row in bed_df.iterrows():
        try:
            center = (row['start'] + row['end']) // 2
            region_start = max(0, center - extend)
            region_end = center + extend
            
            # Extract signals from each replicate group
            group_profiles = []
            for group_handles in all_bw_handles:
                # Average within each replicate group
                replicate_profiles = []
                for bw in group_handles:
                    values = bw.values(row['chr'], region_start, region_end)
                    if values is not None and len(values) > 0:
                        clean_values = np.array([v if v is not None and not np.isnan(v) else 0 for v in values])
                        
                        if len(clean_values) >= n_bins:
                            # Bin by averaging
                            binned_values = []
                            values_per_bin = len(clean_values) // n_bins
                            for i in range(n_bins):
                                start_idx = i * values_per_bin
                                end_idx = (i + 1) * values_per_bin if i < n_bins - 1 else len(clean_values)
                                bin_mean = np.mean(clean_values[start_idx:end_idx])
                                binned_values.append(bin_mean)
                            replicate_profiles.append(binned_values)
                        else:
                            padded_values = np.pad(clean_values, (0, n_bins - len(clean_values)), 'constant')
                            replicate_profiles.append(padded_values[:n_bins])
                    else:
                        replicate_profiles.append([0] * n_bins)
                
                # Average across replicates in this group
                if replicate_profiles:
                    group_consensus = np.mean(replicate_profiles, axis=0)
                    group_profiles.append(group_consensus)
                else:
                    group_profiles.append([0] * n_bins)
            
            # Average across all groups (if multiple conditions)
            if group_profiles:
                final_profile = np.mean(group_profiles, axis=0)
                all_profiles.append(final_profile)
                valid_regions += 1
                
        except Exception as e:
            continue
    
    # Close all bigwig files
    for group_handles in all_bw_handles:
        for bw in group_handles:
            bw.close()
    
    if not all_profiles:
        return None
    
    # Convert to numpy array and calculate mean profile
    all_profiles = np.array(all_profiles)
    mean_profile = np.mean(all_profiles, axis=0)
    std_profile = np.std(all_profiles, axis=0)
    
    # Create position array (relative to center)
    positions = np.linspace(-extend, extend, n_bins)
    
    return {
        'positions': positions,
        'mean_signal': mean_profile,
        'std_signal': std_profile,
        'all_profiles': all_profiles,
        'n_regions': valid_regions
    }

def create_single_boxplot(signals_dict, bigwig_names, bed_names_ordered, y_axis_max, original_bed_names=None):
    """Create boxplot for single or multiple bigwig signals across groups"""
    
    plot_data = []
    
    # Create mapping from custom names to original names if provided
    if original_bed_names is not None:
        name_mapping = dict(zip(bed_names_ordered, original_bed_names))
    else:
        name_mapping = {name: name for name in bed_names_ordered}
    
    # Handle single vs multiple bigwigs
    if len(bigwig_names) == 1:
        # Single bigwig - use ordered bed names
        for custom_name in bed_names_ordered:
            original_name = name_mapping[custom_name]
            if original_name in signals_dict:
                signals = signals_dict[original_name]
                for signal in signals:
                    plot_data.append({
                        'Group': custom_name,  # Use custom name for display
                        'Signal': signal
                    })
        
        df = pd.DataFrame(plot_data)
        available_groups = bed_names_ordered  # Use custom names for ordering
        
        fig, ax = plt.subplots(figsize=(max(8, len(available_groups) * 1.2), 6))
        
        box_plot = sns.boxplot(data=df, x='Group', y='Signal', 
                              order=available_groups, 
                              palette='Set2',
                              showfliers=False,
                              ax=ax)
        
        ax.set_title(f'{bigwig_names[0]} Signal Distribution at Peak Centers\nAcross Uploaded BED Files', 
                    fontsize=14, fontweight='bold', pad=20)
        
    else:
        # Multiple bigwigs - create paired comparison
        for bigwig_idx, bigwig_name in enumerate(bigwig_names):
            signals_dict_for_bigwig = signals_dict[bigwig_idx]
            for custom_name in bed_names_ordered:
                original_name = name_mapping[custom_name]
                if original_name in signals_dict_for_bigwig:
                    signals = signals_dict_for_bigwig[original_name]
                    for signal in signals:
                        plot_data.append({
                            'Group': custom_name,  # Use custom name for display
                            'BigWig': bigwig_name,
                            'Signal': signal
                        })
        
        df = pd.DataFrame(plot_data)
        available_groups = bed_names_ordered  # Use custom names for ordering
        
        # Check if we have data to plot
        if df.empty:
            raise ValueError("No data available for plotting")
        
        fig, ax = plt.subplots(figsize=(max(12, len(available_groups) * 1.5), 6))
        
        # FIXED: Use explicit boxplot parameters to avoid seaborn issues
        try:
            box_plot = sns.boxplot(
                data=df, 
                x='Group', 
                y='Signal', 
                hue='BigWig',
                order=available_groups,
                hue_order=bigwig_names,  # This ensures upload order is preserved!
                palette=['lightblue', 'lightcoral', 'lightgreen', 'lightyellow'][:len(bigwig_names)],
                showfliers=False,
                ax=ax,
                # Add explicit box styling to avoid seaborn internal issues
                boxprops=dict(alpha=0.7),
                whiskerprops=dict(alpha=0.7),
                capprops=dict(alpha=0.7),
                medianprops=dict(alpha=0.7)
            )
        except Exception as e:
            # Fallback: create simpler boxplot without hue
            st.warning(f"Creating simplified boxplot due to: {e}")
            
            # Reshape data for simple plotting
            simple_data = []
            for group_name in available_groups:
                group_data = df[df['Group'] == group_name]
                if not group_data.empty:
                    for bigwig_name in bigwig_names:
                        bigwig_data = group_data[group_data['BigWig'] == bigwig_name]
                        if not bigwig_data.empty:
                            combined_name = f"{group_name}\n({bigwig_name})"
                            for signal in bigwig_data['Signal']:
                                simple_data.append({
                                    'Combined_Group': combined_name,
                                    'Signal': signal
                                })
            
            simple_df = pd.DataFrame(simple_data)
            box_plot = sns.boxplot(
                data=simple_df, 
                x='Combined_Group', 
                y='Signal',
                palette='Set2',
                showfliers=False,
                ax=ax
            )
        
        ax.set_title(f'Signal Distribution Comparison at Peak Centers\nAcross Uploaded BED Files', 
                    fontsize=14, fontweight='bold', pad=20)
        
        # Only add legend if we successfully created the hue plot
        if 'BigWig' in df.columns and len(df['BigWig'].unique()) > 1:
            ax.legend(title='BigWig Groups', loc='upper right')
    
    ax.set_xlabel('BED File Groups', fontsize=12, fontweight='bold')
    ax.set_ylabel('Mean Signal', fontsize=12, fontweight='bold')
    
    # FIX: Use plt.setp instead of ax.set_xticklabels to avoid the warning
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_ylim(0, y_axis_max)
    
    plt.tight_layout()
    return fig

def create_subplot_line_plot(profile_dict_list, bigwig_names, bed_names_ordered, original_bed_names=None):
    """Create line plot with separate subplots for each group - PRESERVING BED FILE ORDER"""
    
    # Create mapping from custom names to original names if provided
    if original_bed_names is not None:
        name_mapping = dict(zip(bed_names_ordered, original_bed_names))
    else:
        name_mapping = {name: name for name in bed_names_ordered}
    
    # FIXED: Check each BED file individually per BigWig group
    valid_groups = []
    for custom_name in bed_names_ordered:
        original_name = name_mapping[custom_name]
        
        # Check if at least one BigWig group has data for this BED file
        has_data_in_any_bigwig = False
        for bigwig_idx in range(len(bigwig_names)):
            if (original_name in profile_dict_list[bigwig_idx] and 
                profile_dict_list[bigwig_idx][original_name] is not None):
                has_data_in_any_bigwig = True
                break
        
        if has_data_in_any_bigwig:
            valid_groups.append((custom_name, original_name))
    
    if not valid_groups:
        st.error("No valid profiles to plot")
        return None
    
    n_groups = len(valid_groups)
    st.write(f"DEBUG: Found {n_groups} valid groups to plot: {[g[0] for g in valid_groups]}")
    
    # Calculate subplot layout
    if n_groups <= 3:
        ncols = n_groups
        nrows = 1
    elif n_groups <= 6:
        ncols = 3
        nrows = 2
    else:
        ncols = 4
        nrows = (n_groups + 3) // 4
    
    # Create figure with subplots - INCREASED HEIGHT FOR BETTER SPACING
    fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 4*nrows + 1.5), sharey=True)
    
    # Handle single subplot case
    if n_groups == 1:
        axes = [axes]
    elif nrows == 1:
        axes = axes if hasattr(axes, '__iter__') else [axes]
    else:
        axes = axes.flatten()
    
    # Colors for different bigwigs - MAINTAIN UPLOAD ORDER
    colors = ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D'][:len(bigwig_names)]
    
    # Plot each group in its own subplot - USING ORDERED GROUPS
    for idx, (custom_name, original_name) in enumerate(valid_groups):
        ax = axes[idx]
        
        st.write(f"DEBUG: Plotting subplot {idx} for {custom_name} (original: {original_name})")
        
        # Plot data for each BigWig that has this BED file
        plotted_any = False
        for bigwig_idx, bigwig_name in enumerate(bigwig_names):
            if (original_name in profile_dict_list[bigwig_idx] and 
                profile_dict_list[bigwig_idx][original_name] is not None):
                
                profile_data = profile_dict_list[bigwig_idx][original_name]
                
                positions = profile_data['positions']
                mean_signal = profile_data['mean_signal']
                std_signal = profile_data['std_signal']
                n_regions = profile_data['n_regions']
                
                st.write(f"DEBUG: Plotting {bigwig_name} for {custom_name}, n_regions={n_regions}")
                
                # Plot mean line - use upload order for line styles
                line_style = '-' if bigwig_idx == 0 else '--'
                ax.plot(positions, mean_signal, 
                        color=colors[bigwig_idx], 
                        linewidth=2,
                        linestyle=line_style,
                        label=f'{bigwig_name} (n={n_regions})')
                
                # Add confidence interval (mean ± SEM)
                sem_signal = std_signal / np.sqrt(n_regions)
                ax.fill_between(positions, 
                               mean_signal - sem_signal, 
                               mean_signal + sem_signal,
                               color=colors[bigwig_idx], 
                               alpha=0.15)
                
                plotted_any = True
        
        if not plotted_any:
            st.warning(f"No data plotted for {custom_name}")
            continue
        
        # Customize individual subplot - use custom name for title
        ax.set_title(f"{custom_name}", fontsize=11, fontweight='bold')
        ax.axvline(x=0, color='black', linestyle=':', alpha=0.5, linewidth=1)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(positions[0], positions[-1])
        
        if len(bigwig_names) > 1:
            ax.legend(fontsize=8, loc='upper right')
        
        # Set x-axis labels
        if idx >= (nrows-1) * ncols:  # Bottom row
            ax.set_xlabel('Distance from Peak Center (bp)', fontsize=10)
        
        # Set y-axis label for leftmost plots
        if idx % ncols == 0:  # Leftmost column
            ax.set_ylabel('Mean Signal', fontsize=10)
    
    # Hide unused subplots
    for idx in range(n_groups, len(axes)):
        axes[idx].set_visible(False)
    
    # Set overall title with better positioning
    if len(bigwig_names) == 1:
        title = f'{bigwig_names[0]} Signal Profiles Across Uploaded BED Files'
    else:
        title = f'Signal Profile Comparison Across Uploaded BED Files\n(Solid=1st Upload, Dashed=2nd Upload, etc.)'
    
    # IMPROVED TITLE POSITIONING
    fig.suptitle(title, fontsize=14, fontweight='bold', y=0.98)
    
    # BETTER SPACING ADJUSTMENT
    plt.tight_layout()
    plt.subplots_adjust(top=0.85 if len(bigwig_names) > 1 else 0.90)
    
    return fig

def main():
    st.title("📊 BigWig Signal Analysis Tool")
    st.markdown("""
    Upload BigWig and BED files to analyze signal distributions.
    - **Boxplots**: Show signal distribution at peak centers (±500bp, customizable)
    - **Line plots**: Show signal profiles across regions (±2000bp, auto-scaled)
    - **Replicate Averaging**: Automatically average biological replicates during signal extraction
    - **Signal Matrix Export/Import**: Export extracted signals and re-import for instant plotting
    """)
    
    # Check if we have pre-extracted data
    pre_extracted_data = None
    
    st.header("🔄 Import Pre-Extracted Signal Data (Optional)")
    st.info("Upload a previously exported Excel file to skip signal extraction and go directly to plotting.")
    
    uploaded_excel = st.file_uploader(
        "Choose exported signal data file:",
        type=['xlsx'],
        help="Upload an Excel file previously exported from this tool",
        key="excel_uploader"
    )
    
    if uploaded_excel:
        with st.spinner("Loading pre-extracted data..."):
            pre_extracted_data = load_signal_data_from_excel(uploaded_excel)
        
        if pre_extracted_data:
            st.success("✅ Successfully loaded pre-extracted data!")
            
            # Show summary
            st.write("**Data Summary:**")
            col1, col2 = st.columns(2)
            
            with col1:
                st.write("**Original BigWig Files:**")
                for i, name in enumerate(pre_extracted_data['bigwig_file_names'], 1):
                    st.write(f"{i}. {name}")
                
                st.write("**BigWig Groups:**")
                for i, name in enumerate(pre_extracted_data['group_names'], 1):
                    st.write(f"{i}. {name}")
            
            with col2:
                st.write("**Original BED Files:**")
                for i, name in enumerate(pre_extracted_data['bed_names_ordered'], 1):
                    st.write(f"{i}. {name}")
                
                st.write("**Analysis Date:**")
                st.write(pre_extracted_data['metadata'].get('Analysis Date', 'Unknown'))
    
    # If we have pre-extracted data, skip to plotting
    if pre_extracted_data:
        st.markdown("---")
        st.header("🎨 Plot Generation from Pre-Extracted Data")
        
        # Custom names for imported data
        custom_bigwig_names, custom_bed_names = setup_custom_names(
            pre_extracted_data['group_names'], 
            pre_extracted_data['bed_names_ordered'],
            mode="imported"
        )
        
        # Settings for plotting
        with st.sidebar:
            st.header("🔧 Plot Settings")
            
            # Use original plot type or allow override
            original_plot_type = pre_extracted_data['analysis_params'].get('plot_type', 'Both')
            plot_type = st.selectbox(
                "Select plot type:",
                ["Boxplot", "Line plot", "Both"],
                index=["Boxplot", "Line plot", "Both"].index(original_plot_type) if original_plot_type in ["Boxplot", "Line plot", "Both"] else 0,
                help="Choose the type of visualization"
            )
            
            # Boxplot settings
            if plot_type in ["Boxplot", "Both"]:
                st.subheader("Boxplot Settings")
                original_y_max = pre_extracted_data['analysis_params'].get('y_max', 25.0)
                y_max = st.number_input(
                    "Y-axis maximum:",
                    min_value=0.1,
                    value=float(original_y_max),
                    step=0.1,
                    help="Maximum value for boxplot y-axis"
                )
            
            # Export format settings
            st.subheader("Export Settings")
            export_format = st.selectbox(
                "Export format:",
                ["PNG", "PDF"],
                help="Choose file format for plot downloads"
            )
        
        signals_data = pre_extracted_data['signals_data']
        profile_data = pre_extracted_data['profile_data']
        
        # Store original names for data lookup
        original_bed_names = pre_extracted_data['bed_names_ordered']
        original_group_names = pre_extracted_data['group_names']
        
        # Store plots in session state to prevent disappearing
        if 'current_plots' not in st.session_state:
            st.session_state.current_plots = {}
        
        # Generate plots instantly
        if plot_type in ["Boxplot", "Both"] and signals_data:
            st.header("📊 Boxplot Results")
            try:
                if len(original_group_names) == 1:
                    signals_dict = signals_data[0]
                else:
                    signals_dict = signals_data
                
                # Create plot with name mapping
                fig_box = create_single_boxplot(
                    signals_dict, 
                    custom_bigwig_names, 
                    custom_bed_names, 
                    y_max, 
                    original_bed_names=original_bed_names
                )
                st.session_state.current_plots['boxplot'] = fig_box
                
                # Display plot
                st.pyplot(fig_box, use_container_width=True)
                
                # Download button with format selection
                plot_data, filename, mime_type = export_plot_with_format(fig_box, "signal_boxplot", export_format)
                
                st.download_button(
                    label=f"📥 Download Boxplot ({export_format})",
                    data=plot_data,
                    file_name=filename,
                    mime=mime_type
                )
                
            except Exception as e:
                st.error(f"Error creating boxplot: {e}")
                st.exception(e)
        
        if plot_type in ["Line plot", "Both"] and profile_data:
            st.header("📈 Line Plot Results")
            try:
                # Create plot with name mapping
                fig_line = create_subplot_line_plot(
                    profile_data, 
                    custom_bigwig_names, 
                    custom_bed_names, 
                    original_bed_names=original_bed_names
                )
                
                if fig_line:
                    st.session_state.current_plots['lineplot'] = fig_line
                    
                    # Display plot
                    st.pyplot(fig_line, use_container_width=True)
                    
                    # Download button with format selection
                    plot_data, filename, mime_type = export_plot_with_format(fig_line, "signal_lineplot", export_format)
                    
                    st.download_button(
                        label=f"📥 Download Line Plot ({export_format})",
                        data=plot_data,
                        file_name=filename,
                        mime=mime_type
                    )
                
            except Exception as e:
                st.error(f"Error creating line plot: {e}")
                st.exception(e)
        
        return  # Skip the rest of the UI
    
    # Original analysis workflow
    st.markdown("---")
    st.header("🔬 New Analysis")
    st.info("Upload BigWig and BED files to perform signal extraction and analysis.")
    
    # Initialize analysis state
    if 'analysis_running' not in st.session_state:
        st.session_state.analysis_running = False
    
    # Sidebar for file uploads and settings
    with st.sidebar:
        st.header("🔧 Settings")
        
        # Plot type selection
        plot_type = st.selectbox(
            "Select plot type:",
            ["Boxplot", "Line plot", "Both"],
            help="Choose the type of visualization",
            disabled=st.session_state.analysis_running  # Disable during analysis
        )
        
        # Boxplot settings
        if plot_type in ["Boxplot", "Both"]:
            st.subheader("Boxplot Settings")
            y_max = st.number_input(
                "Y-axis maximum:",
                min_value=0.1,
                value=25.0,  # Changed default to match your image
                step=0.1,
                help="Maximum value for boxplot y-axis",
                disabled=st.session_state.analysis_running
            )
            
            extend_bp = st.number_input(
                "Signal window (±bp from center):",
                min_value=50,
                max_value=5000,
                value=500,
                step=50,
                help="How many base pairs around peak center to analyze",
                disabled=st.session_state.analysis_running
            )
        
        # File upload settings
        st.subheader("Data Limits")
        max_regions = st.number_input(
            "Maximum regions per BED file:",
            min_value=100,
            max_value=10000,
            value=5000,
            step=100,
            help="Randomly sample if BED file has more regions",
            disabled=st.session_state.analysis_running
        )
        
        # Export format settings
        st.subheader("Export Settings")
        export_format = st.selectbox(
            "Export format:",
            ["PNG", "PDF"],
            help="Choose file format for plot downloads",
            disabled=st.session_state.analysis_running
        )
    
    # Main content area
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.header("📁 Upload BigWig Files")
        bigwig_files = st.file_uploader(
            "Choose BigWig files:",
            type=['bw', 'bigwig'],
            accept_multiple_files=True,
            help="Upload 1-4 BigWig files. Multiple files will be compared.",
            disabled=st.session_state.analysis_running  # Disable during analysis
        )
        
        if bigwig_files:
            st.success(f"Uploaded {len(bigwig_files)} BigWig file(s)")
            for i, bw_file in enumerate(bigwig_files, 1):
                st.write(f"{i}. {bw_file.name}")
            
            # Replicate grouping interface
            replicate_groups = None
            if len(bigwig_files) > 1:
                enable_grouping = st.checkbox(
                    "🔗 Enable replicate grouping",
                    value=False,  # Start with disabled since you have individual files
                    help="Group biological replicates to average their signals",
                    disabled=st.session_state.analysis_running
                )
                
                if enable_grouping:
                    replicate_groups = setup_replicate_groups(bigwig_files)
                else:
                    # Each file is its own group
                    replicate_groups = [[i] for i in range(len(bigwig_files))]
            else:
                replicate_groups = [[0]]  # Single file
    
    with col2:
        st.header("📄 Upload BED Files")
        bed_files = st.file_uploader(
            "Choose BED files:",
            type=['bed'],
            accept_multiple_files=True,
            help="Upload BED files in the order you want them to appear in plots",
            disabled=st.session_state.analysis_running  # Disable during analysis
        )
        
        if bed_files:
            st.success(f"Uploaded {len(bed_files)} BED file(s)")
            for i, bed_file in enumerate(bed_files, 1):
                st.write(f"{i}. {bed_file.name}")
    
    # Analysis button with state management
    if st.session_state.analysis_running:
        st.warning("⏳ Analysis in progress... Please wait.")
        if st.button("❌ Cancel Analysis", type="secondary"):
            st.session_state.analysis_running = False
            st.rerun()
        return
    
    if st.button("🚀 Run Analysis", type="primary", use_container_width=True):
        if not bigwig_files:
            st.error("Please upload at least one BigWig file")
            return
        
        if not bed_files:
            st.error("Please upload at least one BED file")
            return
        
        if len(bigwig_files) > 4:
            st.error("Please upload no more than 4 BigWig files")
            return
        
        if len(bigwig_files) > 1 and replicate_groups is None:
            st.error("Please define replicate groups or disable grouping")
            return
        
        # Set analysis state to prevent reruns
        st.session_state.analysis_running = True
        
        # Create temporary directory for files
        with tempfile.TemporaryDirectory() as temp_dir:
            try:
                # Save uploaded files
                bigwig_paths = []
                bigwig_file_names = []
                for bw_file in bigwig_files:
                    bw_path = save_uploaded_file(bw_file, temp_dir)
                    bigwig_paths.append(bw_path)
                    bigwig_file_names.append(bw_file.name)
                
                bed_paths = []
                bed_names = []
                
                for bed_file in bed_files:
                    bed_path = save_uploaded_file(bed_file, temp_dir)
                    bed_paths.append(bed_path)
                    bed_names.append(Path(bed_file.name).stem)
                
                # Keep bed files in upload order
                bed_names_ordered = bed_names.copy()
                
                # Create grouped file paths for extraction
                grouped_bigwig_paths = []
                group_names = []
                
                for group in replicate_groups:
                    group_paths = [bigwig_paths[i] for i in group]
                    grouped_bigwig_paths.append(group_paths)
                    
                    # Create group name - PRESERVE UPLOAD ORDER
                    if len(group) == 1:
                        group_name = Path(bigwig_files[group[0]].name).stem
                    else:
                        # Use name of first file in group to preserve order
                        first_file_name = Path(bigwig_files[group[0]].name).stem
                        group_name = f"{first_file_name}_group"
                    group_names.append(group_name)
                
                # Progress tracking with containers
                progress_container = st.container()
                with progress_container:
                    progress_bar = st.progress(0)
                    status_text = st.empty()
                
                # Extract signals
                signals_data = None
                if plot_type in ["Boxplot", "Both"]:
                    status_text.text("Extracting signals for boxplots...")
                    
                    signals_data = []
                    total_tasks = len(grouped_bigwig_paths) * len(bed_paths)
                    current_task = 0
                    
                    for group_idx, group_paths in enumerate(grouped_bigwig_paths):
                        signals_dict = {}
                        for i, (bed_path, bed_name) in enumerate(zip(bed_paths, bed_names)):
                            current_task += 1
                            progress = (current_task / total_tasks) * (0.5 if plot_type == "Both" else 1.0)
                            progress_bar.progress(progress)
                            
                            signals = extract_signals_fast([group_paths], bed_path, extend=extend_bp, max_regions=max_regions)
                            if signals:
                                signals_dict[bed_name] = signals
                        signals_data.append(signals_dict)
                
                profile_data = None
                if plot_type in ["Line plot", "Both"]:
                    status_text.text("Extracting profiles for line plots...")
                    
                    profile_data = []
                    total_tasks = len(grouped_bigwig_paths) * len(bed_paths)
                    current_task = 0
                    
                    for group_idx, group_paths in enumerate(grouped_bigwig_paths):
                        profile_dict = {}
                        for i, (bed_path, bed_name) in enumerate(zip(bed_paths, bed_names)):
                            current_task += 1
                            if plot_type == "Both":
                                progress = 0.5 + (current_task / total_tasks) * 0.5
                            else:
                                progress = current_task / total_tasks
                            progress_bar.progress(progress)
                            
                            profile = extract_signals_for_profile([group_paths], bed_path, extend=2000, max_regions=max_regions, bin_size=20)
                            if profile:
                                profile_dict[bed_name] = profile
                        profile_data.append(profile_dict)
                
                progress_bar.progress(1.0)
                status_text.text("Creating plots...")
                
                # Custom names setup
                custom_bigwig_names, custom_bed_names = setup_custom_names(
                    group_names, bed_names_ordered, mode="new_analysis"
                )
                
                # Store analysis parameters for export
                analysis_params = {
                    'plot_type': plot_type,
                    'y_max': y_max if plot_type in ["Boxplot", "Both"] else None,
                    'extend_bp': extend_bp if plot_type in ["Boxplot", "Both"] else None,
                    'max_regions': max_regions,
                    'line_extend': 2000,
                    'line_bin_size': 20
                }
                
                # Store plots in session state to prevent disappearing
                if 'current_plots' not in st.session_state:
                    st.session_state.current_plots = {}
                
                # Create plots
                if plot_type in ["Boxplot", "Both"] and signals_data:
                    st.header("📊 Boxplot Results")
                    try:
                        if len(custom_bigwig_names) == 1:
                            signals_dict = signals_data[0]
                        else:
                            signals_dict = signals_data
                            
                        fig_box = create_single_boxplot(signals_dict, custom_bigwig_names, custom_bed_names, y_max)
                        st.session_state.current_plots['boxplot'] = fig_box
                        
                        st.pyplot(fig_box, use_container_width=True)
                        
                        # Download button with format selection
                        plot_data, filename, mime_type = export_plot_with_format(fig_box, "signal_boxplot", export_format)
                        
                        st.download_button(
                            label=f"📥 Download Boxplot ({export_format})",
                            data=plot_data,
                            file_name=filename,
                            mime=mime_type
                        )
                        
                    except Exception as e:
                        st.error(f"Error creating boxplot: {e}")
                        st.exception(e)
                
                if plot_type in ["Line plot", "Both"] and profile_data:
                    st.header("📈 Line Plot Results")
                    try:
                        fig_line = create_subplot_line_plot(profile_data, custom_bigwig_names, custom_bed_names)
                        
                        if fig_line:
                            st.session_state.current_plots['lineplot'] = fig_line
                            
                            st.pyplot(fig_line, use_container_width=True)
                            
                            # Download button with format selection
                            plot_data, filename, mime_type = export_plot_with_format(fig_line, "signal_lineplot", export_format)
                            
                            st.download_button(
                                label=f"📥 Download Line Plot ({export_format})",
                                data=plot_data,
                                file_name=filename,
                                mime=mime_type
                            )
                        
                    except Exception as e:
                        st.error(f"Error creating line plot: {e}")
                        st.exception(e)
                
                # Export signal data to Excel
                st.header("💾 Export Signal Data")
                st.info("Export extracted signal data to Excel file for future use. This allows you to skip signal extraction and generate plots instantly.")
                
                try:
                    excel_buffer = export_signal_data_to_excel(
                        signals_data, profile_data, custom_bigwig_names, custom_bed_names, 
                        bigwig_file_names, analysis_params
                    )
                    
                    # Generate filename with timestamp
                    timestamp = pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')
                    filename = f"bigwig_signal_data_{timestamp}.xlsx"
                    
                    st.download_button(
                        label="📊 Download Signal Data (Excel)",
                        data=excel_buffer,
                        file_name=filename,
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        help="Download extracted signal data for future analysis"
                    )
                    
                    st.success("✅ Signal data export ready! You can re-upload this file later to skip signal extraction.")
                    
                except Exception as e:
                    st.error(f"Error exporting signal data: {e}")
                    st.exception(e)
                
                status_text.text("✅ Analysis complete!")
                
            except Exception as e:
                st.error(f"An error occurred during analysis: {e}")
                st.exception(e)
            finally:
                # Reset analysis state
                st.session_state.analysis_running = False

if __name__ == "__main__":
    main()