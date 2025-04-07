Namibia Rodent Microbiome Project
This repository contains all scripts, data, and protocols related to the analysis of 16S Nanopore sequencing data from rodents sampled in Namibia. The goal is to investigate host-microbiome interactions across environmental gradients and rodent species.

🧬 Project Overview
Location: Namibian farmland and wildlife settings

Samples: Rodents trapped in 2023

Sequencing: Full-length 16S rRNA gene using Oxford Nanopore

Focus: Taxonomic profiling, diversity analyses, metadata integration

📂 Repository Structure
Namibia_project/
│
├── Data/                     # All project data
│   ├── raw/                 # Raw inputs (field tracking, PCR plates, robot inputs)
│   ├── processed/           # Cleaned & annotated outputs (Nanodrop, EMU results)
│
├── Protocols/               # Laboratory and data analysis documentation
│   ├── Data_processing/     # RMarkdown summaries of taxonomic filtering strategies
│   └── Lab_protocols/       # Wet lab SOPs (e.g., DNA extraction, AMPure cleanup)
│
├── Results/
│   ├── Figures/             # All plots and graphics
│   └── Tables/              # Exported tabular results
│
├── Scripts/
│   ├── Preprocessing/       # Cleaning + metadata integration scripts
│   ├── Visualization/       # Species maps, diversity plots, exploratory figures
│   ├── Basecalling/         # HPC SLURM scripts for Dorado and EMU
│   └── Master_script.R      # Main pipeline coordinator script
│
├── Manuscript/              # Paper drafts (WIP)
├── Meetings/                # Project meeting notes & reports
└── README.md                # You're here!
🔧 Pipeline Summary
The master script (Scripts/Master_script.R) coordinates all core analyses:

🧪 Part 0 – External Preprocessing (HPC)
Basecalling: Dorado

Filtering: Length and quality filtering of 16S reads (NanoPlot summary)

Taxonomy: EMU classifier using SILVA

📋 Part 1–3 – Metadata & QC Integration
Cleaning rodent field metadata

Nanodrop QC and AMPure cleanup assessment

Sample metadata linked to barcode IDs

🧬 Part 4 – OTU + Taxonomy + Filtering Comparison
Merging OTU tables with taxonomy and rodent metadata

Comparing Marly's vs. Melanie's EMU filtering outputs

Alpha diversity, rare taxa, phylum composition, rarefaction, NMDS

Output:
Protocols/Data_processing/04b_compare_filtering_marly_vs_melanie.html
Protocols/Data_processing/Comparison_filtering_strategies_Marly_Melanie.pdf

📊 Part 5 – Visualization
Species maps and scatterplots

Habitat-specific diversity

Sex, age, weight, and temporal patterns

📑 Documentation & Protocols
Wet lab protocols: Protocols/Lab_protocols/

EMU filtering comparison: 04b_compare_filtering_marly_vs_melanie.Rmd

Full documentation of Nanopore preprocessing:
Protocols/Data_processing/README_16s_data_preprocessing.md

👩‍🔬 Authors
Fay Webster 

Marly Erazo

Otto Netzel

Lilla Jordan

Emanuel Heitlinger

Conor Noonan

Dong Xia

Melanie Hay

📌 Getting Started
To run the analysis:

Clone the repo locally

Open Scripts/Master_script.R

Adjust paths if needed via here::here()

Run the script step-by-step or modularly via sourced files


