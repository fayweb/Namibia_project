# Namibia_project
Folders:
Data = storage of data
  - Nanodrop_measurements = Original twk and tsv files 
  - CSV = converted files to csv
  - Rodents_catching_data.xlsx = Excel sheet: Tracking mouse trapping data,
    includes data on coordinates, morphological species identification, sex, weight
    
 R = storage of R code
  - Data_to_input = managing, preparing files from measurements to be R code ready
    - Nanodrop_data_management = Cleaning, merging and analyzing nanodrop data
  - Input_to_product = Code to transform clean data to data products ready for analysis
    - DCM = Code to merge clean nanodrop, PCR and catching/tracking data
      
# Nanopore Sequencing Data Analysis Repository

## Introduction

Welcome to our nanopore sequencing data analysis repository! This repository is part of our project where we collected rodent samples from Namibian farmland in 2023. The goal of this repository is to store and manage various types of sequencing data, including 16S and 18S sequencing, to help us analyze microbiota and other related data. This README will guide you through the repository's structure and provide clear examples, so you can easily upload and organize data, even if you have no prior experience with data management.

## Repository Structure

Our repository is organized into several main folders, each serving a specific purpose. Below is a detailed explanation of the structure and examples of the types of files that go into each folder:

### Main Folders

1. **`data/`**: This folder contains all the data files related to the project. It is subdivided into several subfolders:
   - **`raw_data/`**: Store raw sequencing data here. These are files in their original, unprocessed form. Examples include:
     - `data/raw_data/16S_sequencing_run1.fastq`
     - `data/raw_data/18S_sequencing_sample_A.fastq`
     - `data/raw_data/metadata_run1.csv`
   - **`processed_data/`**: This is where you place data that has been processed or cleaned. Examples include:
     - `data/processed_data/16S_cleaned_sequences.fasta`
     - `data/processed_data/18S_normalized_counts.csv`
   - **`metadata/`**: Metadata files that describe the datasets should be placed here. These files provide essential information about the data, such as sample IDs, experimental conditions, and sequencing parameters. Examples include:
     - `data/metadata/sample_metadata_16S.csv`
     - `data/metadata/run1_metadata.txt`

2. **`code/`**: This folder contains all the code used for data analysis, including scripts and notebooks.
   - **`scripts/`**: Store standalone scripts for processing and analyzing the data here. Examples include:
     - `code/scripts/16S_data_cleaning.py`
     - `code/scripts/18S_taxonomy_assignment.R`
   - **`notebooks/`**: If you use Jupyter or other notebook interfaces, store them here. Notebooks are useful for interactive data exploration. Examples include:
     - `code/notebooks/16S_sequence_analysis.ipynb`
     - `code/notebooks/18S_data_visualization.ipynb`

3. **`results/`**: This folder is for storing the results of your analyses, such as figures, tables, and summary reports. Examples include:
   - **`figures/`**: Store any figures or plots generated during the analysis here.
     - `results/figures/16S_taxonomic_distribution.png`
   - **`tables/`**: Store tables or data summaries here.
     - `results/tables/18S_abundance_table.csv`
   - **`reports/`**: Store any written reports or summaries here.
     - `results/reports/16S_analysis_summary.pdf`

4. **`docs/`**: Any additional documentation related to the project, such as protocols or references, can be placed here. Examples include:
   - `docs/analysis_pipeline.pdf`
   - `docs/sequencing_protocol.docx`

### Example Data Structure

To help you visualize how the data should be organized, here’s an example of what the repository might look like after you’ve added your data and code:


## How to Add Your Data

### Raw Data

1. **16S Sequencing Data**: If you have raw 16S sequencing data files (e.g., `.fastq` files), place them in the `data/raw_data/` folder. Example:
   - `data/raw_data/16S_sequencing_run1.fastq`

2. **18S Sequencing Data**: Similarly, place any raw 18S sequencing data files in the same folder. Example:
   - `data/raw_data/18S_sequencing_sample_A.fastq`

3. **Metadata**: Add any metadata files that describe the samples, sequencing conditions, or other relevant details to the `data/metadata/` folder. Example:
   - `data/metadata/sample_metadata_16S.csv`

### Processed Data

1. **Processed 16S Data**: Once you’ve cleaned or processed the 16S data (e.g., removing low-quality reads), move those files to the `data/processed_data/` folder. Example:
   - `data/processed_data/16S_cleaned_sequences.fasta`

2. **Processed 18S Data**: Similarly, store any processed 18S data here. Example:
   - `data/processed_data/18S_normalized_counts.csv`

### Metadata

1. **Creating Metadata Files**: For each dataset, create a metadata file that includes the following information:
   - **File Name**: Name of the data file.
   - **Sample IDs**: List of sample IDs.
   - **Date Collected**: Date the samples were collected.
   - **Experimental Conditions**: Any relevant experimental conditions (e.g., treatment groups).
   - **Sequencing Details**: Information about the sequencing run (e.g., platform, read length).

   Example Metadata File (`data/metadata/sample_metadata_16S.csv`):


## File Naming Conventions

To keep everything organized, please follow these naming conventions:

- **Descriptive Names**: Use names that reflect the content, such as `16S_sequencing_run1.fastq`.
- **No Spaces**: Avoid spaces in file names; use underscores (`_`) instead.
- **Version Numbers**: If you’re uploading multiple versions of a file, add a version number (e.g., `16S_cleaned_sequences_v1.fasta`).

## Analyzing the Data

Once your data is uploaded, you can start analyzing it using the scripts and notebooks in the `code/` folder. Here are a few examples of the types of analysis you might perform:

1. **16S Sequencing Analysis**: Use the script `16S_data_cleaning.py` to clean the raw sequences and then analyze them with the notebook `16S_sequence_analysis.ipynb`.

2. **18S Sequencing Analysis**: Use the script `18S_taxonomy_assignment.R` to assign taxonomy to the 18S sequences and then visualize the data with `18S_data_visualization.ipynb`.

## How to Edit the README File

As the project evolves, you might need to update this README to reflect new data, code, or analysis steps. Here’s how you can do that:

1. **Edit the README**: Simply click on the README file in the GitHub repository and click the pencil icon to edit it. You can also clone the repository locally, edit the README file using a text editor, and then push the changes back to GitHub.

2. **What to Add**:
- **New Data Descriptions**: If you add new types of data, update the `How to Add Your Data` section with descriptions of where that data should go.
- **New Analysis Steps**: If you develop new analysis scripts or notebooks, update the `Analyzing the Data` section with details on how to use them.
- **New Results**: If you generate new results, update the `results/` section with descriptions of the new outputs.

3. **Version Control**: Always ensure that any significant changes to the README are committed with a descriptive message, so we can track the evolution of the repository.
