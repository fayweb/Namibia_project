# ***********************************************************
# Title: Namibia Rodent Project - Data Processing & Analysis
# Purpose: This script executes all processing & analysis steps
#          for the rodent microbiome sequencing study.
#
# Authors: Fay Webster, Marly Erazo, Otto Netzel, Lilla Jordan,
#          Emanuel Heitlinger, Conor Noonan, Dong Xia, Melanie Hay
#

# ***********************************************************
# ***********************************************************
# ***********************************************************
# Part 0: 16S Nanopore Preprocessing (External HPC Pipeline)
# ***********************************************************
# Note: The raw 16S Nanopore data were preprocessed outside of R
#       using a shell-based pipeline executed on the MPCDF cluster.
#
# 🧪 Step 0.1: Basecalling & Demultiplexing (Dorado)
# --------------------------------------------------
# Tool: ONT Dorado
# Script: Scripts/Basecalling/dorado_basecall_mpcdf.slurm
#
# Key components:
#   - CUDA module loaded via `module load cuda/11.4`
#   - Basecalling with `--min-qscore 10`
#   - Outputs: BAM files and summary stats
#   - Demultiplexing using `--kit-name EXP-PBC096`
#
# 📂 Outputs:
#   - Basecalled reads:         data/raw/ont_fastq/
#   - Summary & BAM:            data/raw/ont_fastq/*.bam, *.txt

# 🧬 Step 0.2: Quality Filtering & Length Trimming
# --------------------------------------------------
# Filtering performed for full-length 16S reads (1400–1800 bp).
# Quality control and stats generated via:
#   - NanoPlot
#   - samtools stats

# Summary metrics (pre- and post-filtering):
#   ▸ Mean read quality:       Q15.6
#   ▸ Number of reads:         ~1.57M
#   ▸ Read length N50:         ~1,393 bp
#   ▸ Percent reads > Q10:     98.1%
#   ▸ Percent reads > Q15:     65.1%
#
# 📊 Quality report:
#   - NanoPlot HTML report (interactive):
#     ▸ data/raw/16S_sequencing/quality/NanoPlot-report.html

# 🧬 Step 0.3: Taxonomic Classification (EMU)
# --------------------------------------------------
# Tool: EMU (Exact Mapping of 16S reads)
# Script: Scripts/Basecalling/emu_script.slurm
#
# Key components:
#   - Conda environment activation (`EMU`)
#   - Looping over demultiplexed FASTQ files
#   - Running `emu abundance` with the SILVA database
#
# 📂 Outputs:
#   - Taxonomic abundance tables: data/processed/ont_emu_abundance/
#   - One `.tsv` file per sample with tax_id counts and relative abundances

# 📄 Documentation:
#   - Full pipeline details, filtering criteria, and HPC command logs:
#     ▸ Protocols/Data_processing/README_16s_data_preprocessing.pdf

# 🧾 Step 0.4: Sample Tracking Metadata (PCR + Pooling)
# --------------------------------------------------
# The following files document the 16S plate setup, robot buffer input,
# and pooling decisions used to prepare samples for sequencing.
# These are **not used in the downstream R analysis**, but archived
# for provenance and reproducibility.

# 📄 Documents:
#   - 16S PCR plate metadata:
#   data/raw/16S_sequencing/16s_PCR/20231123_PCR_final.xlsx
#   - Opentrons robot buffer input:
#         data/raw/16S_sequencing/16s_PCR/pipetingrobot_muster/Buffer
#CSV_19-03-2024_namibia_plate16s.csv
#   - Opentrons robot sample input:
#         data/raw/16S_sequencing/16s_PCR/pipetingrobot_muster/Sample
#CSV_19-03-2024_namibia_plate16s.csv
#   - Sample pooling design:
#                data/raw/16S_sequencing/16s_PCR/pipetingrobot_muster/
#Sample_pooling_CSV_03-04-2024.csv
# ***********************************************************
# Part 1: Set Standard Settings & Load Libraries ----
# ***********************************************************

# Increase max overlaps for better plotting
options(ggrepel.max.overlaps = Inf)

# Set a reproducible seed
set.seed(13102023)

# Load & install required packages using pacman
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(
  tidyverse, janitor, readr, lubridate, ggplot2, phyloseq, vegan,
  corrplot, patchwork, ggrepel, RColorBrewer, pheatmap, caret,
  randomForest, rfUtilities, optimx, ggpubr, FactoMineR, factoextra,
  leaflet, kableExtra, broom, magrittr, data.table, sf, rnaturalearth,
  RColorBrewer, tmap, mapview, cowplot, magick, readxl, qgraph, vegan, ggvenn,
  decontam, microbiome, corncob, scales, ggsci, microViz, ggraph, ggtext,
  ComplexHeatmap, plotly)

#microbiomAnalysis

# ***********************************************************
# Part 2: Define Project File Paths ----
# ***********************************************************

# Dynamically detect the working directory
project_root <- here::here()

# Define primary directories
data_dir       <- file.path(project_root, "data")
raw_data       <- file.path(data_dir, "raw")
  field_raw_tracking <- file.path(raw_data, "Field_tracking_data")
processed_data <- file.path(data_dir, "processed")
metadata_dir   <- file.path(data_dir, "metadata")
results_dir    <- file.path(project_root, "results")
figures_dir    <- file.path(results_dir, "figures")
tables_dir     <- file.path(results_dir, "tables")
scripts_dir    <- file.path(project_root, "scripts")

# ***********************************************************
# Define Column Selection Vectors for Analysis
# ***********************************************************

# 📌 Trapping & Field Metadata Variables (rodent_data)
trapping_vars <- c(
  "Latitude",
  "Longitude",
  "Morphology_species",
  "Sex",
  "Age",
  "Weight_g",
  "Location_type",
  "Date"
)

# ***********************************************************
# 16S Sequencing Variables
# ***********************************************************

# 🧬 16S Barcode & Sequencing Metadata
barcode_vars_16s <- c("barcode_16s", "conc_16s__PCR", "Gene")

# 🧪 16S Taxonomic Assignments
otu_16s_vars <- c(
  "tax_id",
  "count_16s",
  "species",
  "genus",
  "family",
  "order",
  "class",
  "phylum",
  "superkingdom"
)

# 🔬 Combined 16S Selection Vector
microbiome_16s_vars <- c(barcode_vars_16s, otu_16s_vars)

# ***********************************************************
# 28S Sequencing Variables
# ***********************************************************

# 🧬 28S Barcode & Sequencing Metadata
barcode_vars_28s <- c("barcode_28s")
count_vars_28s   <- c("count_28s")
metadata_vars_28s <- c("Date_28s_PCR", "Gene")

# 🧪 28S Taxonomic Assignments (same structure as 16S)
taxonomy_vars_shared <- c("tax_id", "species", "genus",
                          "family", "order", "class",
                          "phylum", "superkingdom")

# 🔬 Combined 28S Selection Vector
otu_28s_vars <- c(barcode_vars_28s, count_vars_28s,
                  taxonomy_vars_shared, metadata_vars_28s)



# ***********************************************************
# Part 3: Data Cleaning - Rodent Field Data ----
# ***********************************************************

#----------------------------------------------------------*
# 3.1: Import & Clean Rodent Field Data
#----------------------------------------------------------*
message("\n🔹 Step 3.1: Cleaning rodent field data...")
if (1) source(file.path(scripts_dir, "preprocessing", "1_import_clean_field_data.R"))

#----------------------------------------------------------*
# 3.1a: Nanodrop DNA Quality Assessment in the Field
#----------------------------------------------------------*
#message("\n🔹 Step 3.1a: Archiving early Nanodrop DNA QC (field)...")
#source(file.path(scripts_dir, "preprocessing", "01a_nanodrop_field_qc.R"))
#----------------------------------------------------------*
# 3.3: AMPure Cleanup QC (Marly)
#----------------------------------------------------------*
#message("\n🔹 Step 3.3: Processing DNA cleanup QC from AMPure protocol...")
#source(file.path(scripts_dir, "preprocessing", "03_dna_cleaning_qc_marly.R"))



# ***********************************************************
# Part 4: Taxonomic Processing - OTU & Phylogenetic Analysis ----
# ***********************************************************
#----------------------------------------------------------*
#----------------------------------------------------------*
# 4.1a: Rename OTU Columns Using Sample Metadata
#----------------------------------------------------------*
message("\n🔹 Step 4.1a: Linking barcode names to rodent data frame")
if (1) source(file.path(scripts_dir, "preprocessing", "02_link_barcode_to_metadata.R"))
#----------------------------------------------------------*
# 4.1b: Renamed EMU OTU Tables - Marly vs Melanie
#----------------------------------------------------------*
# 📄 Documentation:
#   - EMU outputs and filtering description:
#     ▸ Protocols/Data_processing/EMU_outputs_documentation.md
#----------------------------------------------------------*
#----------------------------------------------------------*
# 4.1c: Integrate OTU Tables + Taxonomy + Rodent Metadata
#----------------------------------------------------------*
# Purpose:
#   - Clean and pivot Marly's and Melanie's OTU count tables
#   - Join with EMU taxonomy assignments
#   - Integrate with rodent metadata using barcodes
#   - Save fully annotated OTU tables separately for each filtering strategy
#
# 📄 Script: Scripts/Preprocessing/04a_merge_otu_tables_to_metadata.R
# 📂 Outputs:
#   - otu_taxonomy_metadata_marly.csv    (standard filtering)
#   - otu_taxonomy_metadata_melanie.csv (lenient filtering)
message(
  "\n🔹 Step 4.1c: Integrating OTU counts with taxonomy and rodent metadata...")
if (1) source(
  file.path(scripts_dir, "preprocessing", "04a_merge_otu_tables_to_metadata.R"))
#----------------------------------------------------------*
# 4.1d: Compare EMU Filtering Pipelines
#----------------------------------------------------------*
# This step runs a full RMarkdown comparison between two
# EMU OTU filtering strategies:
#   - Marly's standard filtering
#   - Melanie's lenient filtering
# 📂 Outputs:
#   ▸ HTML report:
# Protocols/Data_processing/04b_compare_filtering_marly_vs_melanie.html
#   ▸ PDF report:
# Protocols/Data_processing/Comparison_filtering_strategies_Marly_Melanie.pdf
#
# 4.1d: Compare filtering strategies and save report to protocols folder
#message("\n🔹
#Step 4.1d: Comparing filtering approaches from Marly and Melanie...")
#rmarkdown::render(
 # input = file.path(scripts_dir, "preprocessing",
#"04b_compare_filtering_marly_vs_melanie.Rmd"),
  #knit_root_dir = project_root,
#  output_file = "04b_compare_filtering_marly_vs_melanie.html",
 # output_dir = file.path(project_root, "Protocols", "Data_processing")
#)
#----------------------------------------------------------*
# 4.1e: Select Final Filtering Strategy for Downstream Analysis
#----------------------------------------------------------*
# Based on comparisons in 04b_compare_filtering_marly_vs_melanie.Rmd,
# we selected Marly's standard filtering for all downstream microbiome analyses.
# This script copies the final annotated OTU table and standardizes
# its name for clarity and reproducibility.
# 📄 Script: Scripts/Preprocessing/04c_select_final_filtering_output.R
# 📂 Outputs:
#   - otu_taxonomy_metadata_final.csv
#     (used for diversity, phyloseq, and downstream community analysis)
message("\n🔹 Step 4.1e: Finalizing selection of
        Marly's filtered EMU output...")
if (1) source(file.path(scripts_dir, "preprocessing",
                 "04c_select_final_filtering_output.R"))
#----------------------------------------------------------*
# 4.2: Integrate EMU 28S Output from Emanuel Heitlinger
#----------------------------------------------------------*
# Purpose:
#   - Load 28S EMU output and sample metadata from Emanuel
#   - Filter for rodent samples only
#   - Join with rodent metadata
#   - Annotate and prepare for parasite-focused analysis
#
# 📄 Script: Scripts/Preprocessing/05a_mouse_parasite_28s_analysis.R
# 📂 Outputs:
#   - rodent_data with merged 28S data
#   - Used later for parasite diversity and re-clustering
message("\n🔹 Step 4.2: Annotating 28S EMU output from Emanuel
        and joining to rodent metadata...")
if (1) source(file.path(scripts_dir, "preprocessing",
                 "05a_mouse_parasite_28s_analysis.R"))


# ***********************************************************
# Part 5: Microbiome Community Analysis - 16S & 28S ----
# ***********************************************************

# 5.1: Construct Phyloseq Object (16S) ------------------------
if (0) source(file.path(scripts_dir, "Analysis", "0_construct_phyloseq.R"))

# 5.2: Decontamination (Blanks & Controls) -------------------
if (0) source(file.path(scripts_dir, "Analysis", "1_decontamination_separates.R"))

# 5.3: Preprocessing / Filtering -----------------------------
if (0) source(file.path(scripts_dir, "Analysis", "2_quality_summary_stats.R"))

# 5.4: Alpha Diversity ----------------------------------------
if (1) source(file.path(scripts_dir, "Analysis", "3_alpha_diversity.R"))

# 5.5: Beta Diversity -----------------------------------------
if (0) source(file.path(scripts_dir, "Analysis", "04a_beta_distance_ordination.R"))
source(file.path(scripts_dir, "Analysis", "04b_beta_pca_clr_taxa.R"))

# 5.6: Taxonomic Barplots -------------------------------------
if (0) source(file.path(scripts_dir, "Analysis", "5_compositional_plots.R"))

# 5.7: Heatmap of Taxa ----------------------------------------
if (0) source(file.path(scripts_dir, "Analysis", "6_heatmap_microbiome.R"))


































# ***********************************************************
# Part 5: Generate Figures & Reports ----
# ***********************************************************
#----------------------------------------------------------*
# 5.1: Generate Species Distribution Plots
#----------------------------------------------------------*
message("\n🔹 Step 5.1: Generating species distribution plots...")
source(file.path(scripts_dir, "visualization", "05_species_distribution.R"))

#----------------------------------------------------------*
# 5.2: Generate Geographical Mapping
#----------------------------------------------------------*
message("\n🔹 Step 5.2: Mapping species distributions...")
source(file.path(scripts_dir, "visualization", "05_geographical_mapping.R"))

#----------------------------------------------------------*
# 5.3: Analyze Sex & Age Distribution
#----------------------------------------------------------*
message("\n🔹 Step 5.3: Analyzing sex & age distribution...")
source(file.path(scripts_dir, "visualization", "05_sex_age_analysis.R"))
# ----------------------------------------------------------*
# 5.4: Habitat Distribution Analysis
# ----------------------------------------------------------*
message("\n🔹 Step 5.4: Analyzing species distribution by habitat...")
source(file.path(scripts_dir, "visualization", "05_habitat_distribution.R"))

# ----------------------------------------------------------*
# 5.5: Temporal Patterns Analysis
# ----------------------------------------------------------*
message("\n🔹 Step 5.5: Analyzing temporal patterns of species occurrences...")
source(file.path(scripts_dir, "visualization", "05_temporal_patterns.R"))

# ----------------------------------------------------------*
# 5.6: Weight Distribution Analysis
# ----------------------------------------------------------*
message("\n🔹 Step 5.6: Analyzing weight distribution across species...")
source(file.path(scripts_dir, "visualization", "05_weight_distribution.R"))

# ----------------------------------------------------------*
# 5.7: Variable Relationships (Geographical Scatter Plot)
# ----------------------------------------------------------*
message("\n🔹 Step 5.7: Analyzing relationships between species and geographical coordinates...")
source(file.path(scripts_dir, "visualization", "05_variable_relationships.R"))

# ----------------------------------------------------------*
# 5.8: Site Panel Figures (Map Combos)
# ----------------------------------------------------------*
message("\n🔹 Step 5.8: Creating panel plots for site locations...")
source(file.path(scripts_dir, "visualization", "05_site_panels.R"))
# ----------------------------------------------------------*
# 5.9: Mouse Parasite Heatmap - 28S EMU (Emanuel's Analysis)
# ----------------------------------------------------------*
message("\n🔹 Step 5.9: Generating 28S parasite heatmap for rodent samples (Emanuel's EMU output)...")
source(file.path(scripts_dir, "preprocessing", "05a_mouse_parasite_28s_analysis.R"))




# ***********************************************************
# Part 6: Save Final Processed Data ----
# ***********************************************************

#----------------------------------------------------------*
# 6.1: Save Final Processed Data for Further Use
#----------------------------------------------------------*
#message("\n🔹 Step 6.1: Saving final processed data...")
#source(file.path(scripts_dir, "output", "06_save_final_data.R"))

# ***********************************************************
# Part 7: Manuscript Preparation ----
# ***********************************************************

#----------------------------------------------------------*
# 7.1: Compile Report and Manuscript
#----------------------------------------------------------*
#message("\n🔹 Step 7.1: Compiling manuscript and final report...")
#source(file.path(scripts_dir, "manuscript", "07_compile_manuscript.R"))

# ***********************************************************
# Completion Message
# ***********************************************************
message("\n✅ Namibia Rodent Project: All processing and analysis steps completed successfully! 🚀")

