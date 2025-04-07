#----------------------------------------------------------*
# 4.1b: Renamed EMU OTU Tables - Marly vs Melanie
#----------------------------------------------------------*
# Description:
#   Two sets of EMU taxonomic classification outputs are available:
#   - Marly's standard filtering parameters
#   - Melanie's lenient filtering parameters
#
# For each, the following files are provided:
#   - otu_counts_<method>.tsv: Taxonomic abundance (count per tax_id)
#   - taxonomy_<method>.tsv: Taxonomic assignments per tax_id
#
# ✅ These have already been renamed using the cleaned barcode-to-sample mapping
#    and are ready for downstream integration with metadata and use in Phyloseq.
#
# File structure:
#   ▸ data/processed/EMU_output/marly_standard_filtering/
#       ├── otu_counts_marly_standard.tsv
#       └── taxonomy_marly_standard.tsv
#
#   ▸ data/processed/EMU_output/melanie_lenient_filtering/
#       ├── otu_counts_melanie_lenient.tsv
#       └── taxonomy_melanie_lenient.tsv
#
# No additional renaming is required at this step.
# Comparison and selection between filtering strategies
# can be performed in subsequent analysis (e.g., alpha/beta diversity).
#----------------------------------------------------------*
