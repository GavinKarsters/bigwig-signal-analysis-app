# BigWig Signal Analysis Tool

A Streamlit web application for analyzing ChIP-seq signal data from BigWig files across genomic regions defined in BED files.

## Live App

**Access the tool here:** [BigWig Signal Analysis App](https://bigwig-signal-analysis-app-g3gmjpyeepfrz7wrd8pwcn.streamlit.app/)

## Features

### Multiple Visualization Types
- **Boxplots**: Signal distribution analysis at peak centers (customizable window size)
- **Line Plots**: Signal profiles across genomic regions (±2000bp around centers)
- **Comparative Analysis**: Side-by-side comparison of multiple BigWig files

### Flexible Input Options
- **BigWig Files**: Upload 1-4 BigWig files for single or comparative analysis
- **BED Files**: Upload multiple BED files in your desired display order
- **File Format Support**: `.bw`, `.bigwig`, `.bed` formats

### Signal Matrix Export/Import
- **Export Data**: Save extracted signal matrices to Excel files for future use
- **Instant Plotting**: Re-upload exported files to skip signal extraction and generate plots instantly

### Replicate Averaging
- **Group Replicates**: Combine biological replicates during signal extraction
- **Flexible Grouping**: Manual or automatic grouping options
- **Consensus Signals**: Average replicate signals for cleaner visualizations

### Customizable Settings
- **Signal Window**: Adjust base pair window around peak centers 
- **Y-axis Control**: Set maximum values for consistent boxplot scaling
- **Region Sampling**: Control maximum regions per BED file (100-10,000)
- **High Resolution**: Download plots at 300 DPI for publication quality

## How to Use

### Option 1: New Analysis

#### Step 1: Upload Files
1. **BigWig Files**: Upload your signal files (ChIP-seq, ATAC-seq, etc.)
   - Single file: Individual analysis
   - Multiple files: Comparative analysis with different colors/styles
2. **BED Files**: Upload genomic regions of interest
   - Order matters: Files appear in plots in upload order
   - Supports any genomic coordinates in standard BED format

#### Step 2: Configure Settings
- **Plot Type**: Choose boxplot, line plot, or both
- **Replicate Grouping**: Enable to average biological replicates
- **Y-axis Maximum**: Set consistent scaling for boxplots
- **Signal Window**: Adjust analysis window around peak centers
- **Max Regions**: Control computational load and sampling

#### Step 3: Run Analysis
- Click "Run Analysis"
- Monitor progress with real-time status updates
- View interactive results directly in browser

#### Step 4: Export & Download
- **Download Plots**: High-resolution PNG downloads at 300 DPI
- **Export Signal Data**: Save extracted signals to Excel for future use

### Option 2: Import Pre-Extracted Data

#### Step 1: Upload Excel File
- Upload previously exported signal data file
- Automatically loads original analysis parameters and file information

#### Step 2: Instant Plotting
- Adjust plot settings in real-time
- Generate plots instantly without signal extraction
- Download publication-ready figures immediately
