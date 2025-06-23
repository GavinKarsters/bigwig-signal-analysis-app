# BigWig Signal Analysis Tool

An interactive Streamlit web application for analyzing ChIP-seq signal data from BigWig files across genomic regions defined in BED files.

## Live App

**Access the tool here:** [BigWig Signal Analysis App](https://bigwig-signal-analysis-app-g3gmjpyeepfrz7wrd8pwcn.streamlit.app/)

## Features

### Multiple Visualization Types
- **Boxplots**: Signal distribution analysis at peak centers (customizable window size)
- **Line Plots**: Signal profiles across genomic regions (±2000bp around centers)
- **Heatmaps**: Visual comparison of signal intensity across all individual genomic regions, sorted by signal strength.

### Replicate Averaging
- **Group Replicates**: Combine biological replicates during signal extraction
- **Flexible Grouping**: Manual or automatic grouping options
- **Consensus Signals**: Average replicate signals for cleaner visualizations

### Signal Matrix Export/Import
- **Export Data**: Save extracted signal matrices to Excel files for future use
- **Instant Plotting**: Re-upload exported files to skip signal extraction and generate plots instantly



## Workflows


There are two main ways to use the app:

### Option 1: Run a New Analysis
This is the standard workflow for processing raw BigWig and BED files.

1.  **Upload Files**: In the main panel, upload your BigWig and BED files. The upload order determines the plotting order.
2.  **Configure Settings**: Use the sidebar to select the desired plot types, configure replicate grouping, and adjust plot-specific parameters (e.g., heatmap color scale, boxplot y-axis).
3.  **Run Analysis**: Click the `Run Analysis` button.
4.  **View & Export**: The app will display the results. You can then download individual plots or the complete signal data as an Excel file for future use.

### Option 2: Re-Plot from Exported Data
If you have already run an analysis and saved the Excel output, you can use this workflow to save time.

1.  **Upload Excel File**: In the "Import Pre-Extracted Signal Data" section, upload your previously exported `.xlsx` file.
2.  **View & Customize**: The app will instantly generate all plots based on the (meta)data in the file.
3.  **Adjust & Download**: You can use the sidebar to change visualization settings (e.g., switch the heatmap colormap, adjust axis limits) and download the newly styled plots.


## Example output plots

### Boxplot

![Boxplot Example](https://github.com/user-attachments/assets/2df4f3a0-6eb0-4016-ae4d-d2e68adf01de)

### Lineplot

![Lineplot Example](https://github.com/user-attachments/assets/64967602-a0df-4e1f-8645-878ffde48af9)

### Heatmap

![Lineplot Example](https://github.com/user-attachments/assets/11f69cb7-5bb9-4eb0-b689-71fb9a71d628)