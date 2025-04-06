# ***********************************************************
# Title: Namibia Rodent Project - Data Processing & Analysis
# Purpose: This script executes all processing & analysis steps
#          for the rodent microbiome sequencing study.
#
# Authors: Fay Webster, Marly Erazo, Otto Netzel, Lilla Jordan,
#          Emanuel Heitlinger, Conor Noonan, Dong Xia, Melanie Hay
#
# Date: 2024-10-15
# ***********************************************************
# ***********************************************************
# Note: 16S Nanopore Sequencing Preprocessing (Shell-based)
#
# The raw 16S MinION data were preprocessed outside of R using
# a pipeline built from ONT's Dorado and EMU tools on the HPC.
#
# Steps:
#   - Basecalling (Dorado)
#   - Length filtering (1400–1800 bp)
#   - Demultiplexing (Dorado demux)
#   - Quality control (NanoPlot, samtools stats)
#   - Read count summaries
#   - Taxonomic assignment (EMU)
#
# Full protocols and bash commands are described here:
#   ▸ Namibia_project/Protocols/Data_processing/README_16s_data_preprocessing.pdf
#
# Outputs used in R:
#   ▸ data/raw/ont_fastq/
#   ▸ data/processed/ont_emu_abundance/

# ***********************************************************
# Additional Slurm script used for Dorado basecalling on MPCDF:
#   ▸ Scripts/Basecalling/dorado_basecall_mpcdf.slurm
#
# This batch script includes:
#   - Module loading (CUDA)
#   - Dorado basecalling with `--min-qscore 10`
#   - Generation of BAM and summary stats
#   - Demultiplexing with EXP-PBC096 kit

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
  RColorBrewer, tmap, mapview, cowplot, magick, readxl
)

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
# Part 3: Data Cleaning - Rodent Field Data ----
# ***********************************************************

#----------------------------------------------------------*
# 3.1: Import & Clean Rodent Field Data
#----------------------------------------------------------*
message("\n🔹 Step 3.1: Cleaning rodent field data...")
source(file.path(scripts_dir, "preprocessing", "1_import_clean_field_data.R"))

#----------------------------------------------------------*
# 3.1a: Nanodrop DNA Quality Assessment in the Field
#----------------------------------------------------------*
#message("\n🔹 Step 3.1a: Archiving early Nanodrop DNA QC (field)...")
#source(file.path(scripts_dir, "preprocessing", "01a_nanodrop_field_qc.R"))
#----------------------------------------------------------*
# 3.3: AMPure Cleanup QC (Marly)
#----------------------------------------------------------*
message("\n🔹 Step 3.3: Processing DNA cleanup QC from AMPure protocol...")
source(file.path(scripts_dir, "preprocessing", "03_dna_cleaning_qc_marly.R"))



# ***********************************************************
# Part 4: Taxonomic Processing - OTU & Phylogenetic Analysis ----
# ***********************************************************

#----------------------------------------------------------*
# 4.1: Process OTU Table & Taxonomic Assignments
#----------------------------------------------------------*
#message("\n🔹 Step 4.1: Processing OTU & taxonomic data...")
#source(file.path(scripts_dir, "taxonomy", "03_process_OTU_taxonomy.R"))

#----------------------------------------------------------*
# 4.2: Perform Diversity Analysis
#----------------------------------------------------------*
#message("\n🔹 Step 4.2: Running diversity analysis...")
#source(file.path(scripts_dir, "diversity", "04_diversity_analysis.R"))

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

