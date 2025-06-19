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

def extract_signals_fast(bigwig_files, bed_file, extend=500, max_regions=5000):
    """Extract signals for boxplots - peak center only"""
    
    if isinstance(bigwig_files, str):
        bigwig_files = [bigwig_files]
    
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
    
    # Open bigwig files
    try:
        bw_handles = [pyBigWig.open(bw_file) for bw_file in bigwig_files]
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
            
            replicate_signals = []
            for bw in bw_handles:
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
            
            if replicate_signals:
                consensus_signal = np.mean(replicate_signals)
                signals.append(consensus_signal)
                valid_regions += 1
            else:
                signals.append(0)
                
        except Exception as e:
            signals.append(0)
    
    # Close bigwig files
    for bw in bw_handles:
        bw.close()
        
    return signals

def extract_signals_for_profile(bigwig_files, bed_file, extend=2000, max_regions=5000, bin_size=20):
    """Extract signals for line plots - returns position-wise signals"""
    
    if isinstance(bigwig_files, str):
        bigwig_files = [bigwig_files]
    
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
    
    # Open bigwig files
    try:
        bw_handles = [pyBigWig.open(bw_file) for bw_file in bigwig_files]
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
            
            # Extract signals from all replicates
            replicate_profiles = []
            for bw in bw_handles:
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
            
            # Average across replicates
            if replicate_profiles:
                consensus_profile = np.mean(replicate_profiles, axis=0)
                all_profiles.append(consensus_profile)
                valid_regions += 1
                
        except Exception as e:
            continue
    
    # Close bigwig files
    for bw in bw_handles:
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

def create_single_boxplot(signals_dict, bigwig_names, y_axis_max):
    """Create boxplot for single or multiple bigwig signals across groups"""
    
    plot_data = []
    
    # Handle single vs multiple bigwigs
    if len(bigwig_names) == 1:
        # Single bigwig
        for group_name, signals in signals_dict.items():
            for signal in signals:
                plot_data.append({
                    'Group': group_name,
                    'Signal': signal
                })
        
        df = pd.DataFrame(plot_data)
        available_groups = list(signals_dict.keys())
        
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
            for group_name, signals in signals_dict[bigwig_idx].items():
                for signal in signals:
                    plot_data.append({
                        'Group': group_name,
                        'BigWig': bigwig_name,
                        'Signal': signal
                    })
        
        df = pd.DataFrame(plot_data)
        available_groups = list(signals_dict[0].keys())
        
        fig, ax = plt.subplots(figsize=(max(12, len(available_groups) * 1.5), 6))
        
        box_plot = sns.boxplot(data=df, x='Group', y='Signal', hue='BigWig',
                              order=available_groups, 
                              palette=['lightblue', 'lightcoral', 'lightgreen', 'lightyellow'][:len(bigwig_names)],
                              showfliers=False,
                              ax=ax)
        
        ax.set_title(f'Signal Distribution Comparison at Peak Centers\nAcross Uploaded BED Files', 
                    fontsize=14, fontweight='bold', pad=20)
        ax.legend(title='BigWig Files', loc='upper right')
    
    ax.set_xlabel('BED File Groups', fontsize=12, fontweight='bold')
    ax.set_ylabel('Mean Signal', fontsize=12, fontweight='bold')
    ax.set_xticklabels(available_groups, rotation=45, ha='right')
    ax.set_ylim(0, y_axis_max)
    
    plt.tight_layout()
    return fig

def create_subplot_line_plot(profile_dict_list, bigwig_names):
    """Create line plot with separate subplots for each group"""
    
    # Get group names from first bigwig (assuming all have same groups)
    if len(bigwig_names) == 1:
        valid_profiles = {k: v for k, v in profile_dict_list[0].items() if v is not None}
    else:
        # Get groups that exist in all bigwigs
        all_groups = set(profile_dict_list[0].keys())
        for profile_dict in profile_dict_list[1:]:
            all_groups = all_groups.intersection(set(profile_dict.keys()))
        valid_profiles = {group: profile_dict_list[0][group] for group in all_groups 
                         if all(profile_dict_list[i][group] is not None for i in range(len(profile_dict_list)))}
    
    if not valid_profiles:
        st.error("No valid profiles to plot")
        return None
    
    n_groups = len(valid_profiles)
    
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
    
    # Create figure with subplots
    fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 4*nrows + 1), sharey=True)
    
    # Handle single subplot case
    if n_groups == 1:
        axes = [axes]
    elif nrows == 1:
        axes = axes if hasattr(axes, '__iter__') else [axes]
    else:
        axes = axes.flatten()
    
    # Colors for different bigwigs
    colors = ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D'][:len(bigwig_names)]
    
    # Plot each group in its own subplot
    for idx, group_name in enumerate(valid_profiles.keys()):
        ax = axes[idx]
        
        for bigwig_idx, bigwig_name in enumerate(bigwig_names):
            profile_data = profile_dict_list[bigwig_idx][group_name]
            
            if profile_data is not None:
                positions = profile_data['positions']
                mean_signal = profile_data['mean_signal']
                std_signal = profile_data['std_signal']
                n_regions = profile_data['n_regions']
                
                # Plot mean line
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
        
        # Customize individual subplot
        ax.set_title(f"{group_name}", fontsize=11, fontweight='bold')
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
    
    # Set overall title
    if len(bigwig_names) == 1:
        title = f'{bigwig_names[0]} Signal Profiles Across Uploaded BED Files'
    else:
        title = f'Signal Profile Comparison Across Uploaded BED Files\n(Solid=1st BigWig, Dashed=2nd BigWig, etc.)'
    
    fig.suptitle(title, fontsize=16, fontweight='bold', y=0.95)
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.88)
    
    return fig

def main():
    st.title("📊 BigWig Signal Analysis Tool")
    st.markdown("""
    Upload BigWig and BED files to analyze signal distributions.
    - **Boxplots**: Show signal distribution at peak centers (±500bp, customizable)
    - **Line plots**: Show signal profiles across regions (±2000bp, auto-scaled)
    """)
    
    # Sidebar for file uploads and settings
    with st.sidebar:
        st.header("🔧 Settings")
        
        # Plot type selection
        plot_type = st.selectbox(
            "Select plot type:",
            ["Boxplot", "Line plot", "Both"],
            help="Choose the type of visualization"
        )
        
        # Boxplot settings
        if plot_type in ["Boxplot", "Both"]:
            st.subheader("Boxplot Settings")
            y_max = st.number_input(
                "Y-axis maximum:",
                min_value=0.1,
                value=100.0,
                step=0.1,
                help="Maximum value for boxplot y-axis"
            )
            
            extend_bp = st.number_input(
                "Signal window (±bp from center):",
                min_value=50,
                max_value=5000,
                value=500,
                step=50,
                help="How many base pairs around peak center to analyze"
            )
        
        # File upload settings
        st.subheader("Data Limits")
        max_regions = st.number_input(
            "Maximum regions per BED file:",
            min_value=100,
            max_value=10000,
            value=5000,
            step=100,
            help="Randomly sample if BED file has more regions"
        )
    
    # Main content area
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.header("📁 Upload BigWig Files")
        bigwig_files = st.file_uploader(
            "Choose BigWig files:",
            type=['bw', 'bigwig'],
            accept_multiple_files=True,
            help="Upload 1-4 BigWig files. Multiple files will be compared."
        )
        
        if bigwig_files:
            st.success(f"Uploaded {len(bigwig_files)} BigWig file(s)")
            for i, bw_file in enumerate(bigwig_files, 1):
                st.write(f"{i}. {bw_file.name}")
    
    with col2:
        st.header("📄 Upload BED Files")
        bed_files = st.file_uploader(
            "Choose BED files:",
            type=['bed'],
            accept_multiple_files=True,
            help="Upload BED files in the order you want them to appear in plots"
        )
        
        if bed_files:
            st.success(f"Uploaded {len(bed_files)} BED file(s)")
            for i, bed_file in enumerate(bed_files, 1):
                st.write(f"{i}. {bed_file.name}")
    
    # Analysis button
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
        
        # Create temporary directory for files
        with tempfile.TemporaryDirectory() as temp_dir:
            try:
                # Save uploaded files
                bigwig_paths = []
                bigwig_names = []
                
                for bw_file in bigwig_files:
                    bw_path = save_uploaded_file(bw_file, temp_dir)
                    bigwig_paths.append(bw_path)
                    bigwig_names.append(Path(bw_file.name).stem)
                
                bed_paths = []
                bed_names = []
                
                for bed_file in bed_files:
                    bed_path = save_uploaded_file(bed_file, temp_dir)
                    bed_paths.append(bed_path)
                    bed_names.append(Path(bed_file.name).stem)
                
                # Progress bar
                progress_bar = st.progress(0)
                status_text = st.empty()
                
                # Extract signals
                if plot_type in ["Boxplot", "Both"]:
                    status_text.text("Extracting signals for boxplots...")
                    
                    if len(bigwig_paths) == 1:
                        # Single bigwig analysis
                        signals_dict = {}
                        for i, (bed_path, bed_name) in enumerate(zip(bed_paths, bed_names)):
                            progress_bar.progress((i + 1) / len(bed_paths) * 0.5)
                            signals = extract_signals_fast(bigwig_paths, bed_path, extend=extend_bp, max_regions=max_regions)
                            if signals:
                                signals_dict[bed_name] = signals
                        
                        signals_data = signals_dict
                    else:
                        # Multiple bigwig analysis
                        signals_data = []
                        for bw_idx, bw_path in enumerate(bigwig_paths):
                            signals_dict = {}
                            for i, (bed_path, bed_name) in enumerate(zip(bed_paths, bed_names)):
                                progress_bar.progress(((bw_idx * len(bed_paths)) + i + 1) / (len(bigwig_paths) * len(bed_paths)) * 0.5)
                                signals = extract_signals_fast([bw_path], bed_path, extend=extend_bp, max_regions=max_regions)
                                if signals:
                                    signals_dict[bed_name] = signals
                            signals_data.append(signals_dict)
                
                if plot_type in ["Line plot", "Both"]:
                    status_text.text("Extracting profiles for line plots...")
                    
                    profile_data = []
                    for bw_idx, bw_path in enumerate(bigwig_paths):
                        profile_dict = {}
                        for i, (bed_path, bed_name) in enumerate(zip(bed_paths, bed_names)):
                            if plot_type == "Both":
                                progress_bar.progress(0.5 + ((bw_idx * len(bed_paths)) + i + 1) / (len(bigwig_paths) * len(bed_paths)) * 0.5)
                            else:
                                progress_bar.progress(((bw_idx * len(bed_paths)) + i + 1) / (len(bigwig_paths) * len(bed_paths)))
                            
                            profile = extract_signals_for_profile([bw_path], bed_path, extend=2000, max_regions=max_regions, bin_size=20)
                            if profile:
                                profile_dict[bed_name] = profile
                        profile_data.append(profile_dict)
                
                progress_bar.progress(1.0)
                status_text.text("Creating plots...")
                
                # Create plots
                if plot_type in ["Boxplot", "Both"]:
                    st.header("📊 Boxplot Results")
                    try:
                        if len(bigwig_paths) == 1:
                            fig_box = create_single_boxplot(signals_data, bigwig_names, y_max)
                        else:
                            fig_box = create_single_boxplot(signals_data, bigwig_names, y_max)
                        
                        st.pyplot(fig_box)
                        plt.close(fig_box)
                        
                        # Download button for boxplot
                        buf = io.BytesIO()
                        fig_box = create_single_boxplot(signals_data, bigwig_names, y_max) if len(bigwig_paths) == 1 else create_single_boxplot(signals_data, bigwig_names, y_max)
                        fig_box.savefig(buf, format='png', dpi=300, bbox_inches='tight')
                        buf.seek(0)
                        plt.close(fig_box)
                        
                        st.download_button(
                            label="📥 Download Boxplot",
                            data=buf,
                            file_name="signal_boxplot.png",
                            mime="image/png"
                        )
                        
                    except Exception as e:
                        st.error(f"Error creating boxplot: {e}")
                
                if plot_type in ["Line plot", "Both"]:
                    st.header("📈 Line Plot Results")
                    try:
                        fig_line = create_subplot_line_plot(profile_data, bigwig_names)
                        
                        if fig_line:
                            st.pyplot(fig_line)
                            plt.close(fig_line)
                            
                            # Download button for line plot
                            buf = io.BytesIO()
                            fig_line = create_subplot_line_plot(profile_data, bigwig_names)
                            fig_line.savefig(buf, format='png', dpi=300, bbox_inches='tight')
                            buf.seek(0)
                            plt.close(fig_line)
                            
                            st.download_button(
                                label="📥 Download Line Plot",
                                data=buf,
                                file_name="signal_lineplot.png",
                                mime="image/png"
                            )
                        
                    except Exception as e:
                        st.error(f"Error creating line plot: {e}")
                
                status_text.text("✅ Analysis complete!")
                
            except Exception as e:
                st.error(f"An error occurred during analysis: {e}")
                st.exception(e)

if __name__ == "__main__":
    main()