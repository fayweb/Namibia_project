# ***********************************************************
# Title: Merge OTU Counts from EMU to Rodent Metadata
# Purpose: Integrate Marly's and Melanie's filtered OTU tables into rodent_data
# Author: Fay Webster
# ***********************************************************
# This script joins OTU counts from Marly's and Melanie's EMU outputs
# into a single long-format table for easy comparison.
# Includes:
#   - tax_id (unique taxa)
#   - barcode (sequencing ID)
#   - count_marly (standard filtering)
#   - count_melanie (lenient filtering)
# ----------------------------------------------------------*
# Inputs:
#   - OTU count tables:     data/processed/EMU_output/[...]/otu_counts_*.tsv
#   - Taxonomy tables:      data/processed/EMU_output/[...]/taxonomy_*.tsv
#   - Rodent metadata (already integrated with barcode info)
# ----------------------------------------------------------*

message("📥 Loading EMU OTU tables from Marly and Melanie...")

# Define paths
otu_marly_path <- file.path(processed_data, "EMU_output", "marly_standard_filtering", "otu_counts_marly_standard.tsv")
otu_melanie_path <- file.path(processed_data, "EMU_output", "melanie_lenient_filtering", "otu_counts_melanie_lenient.tsv")

# Read tables
otu_marly   <- read_tsv(otu_marly_path, show_col_types = FALSE)
otu_melanie <- read_tsv(otu_melanie_path, show_col_types = FALSE)

# ✅ Keep ALL barcodes (including controls and replicates)
valid_barcodes <- rodent_data %>%
  filter(!is.na(barcode)) %>%
  distinct(barcode) %>%
  pull()

# 📊 Pivot OTU tables to long format and match barcodes
otu_marly_long <- otu_marly %>%
  rename_with(~ str_remove(.x, "e17a8f2887894f8d7becdbeaafbc97db14bc8e66_EXP-PBC096_")) %>%
  pivot_longer(-tax_id, names_to = "barcode", values_to = "count_marly") %>%
  filter(barcode %in% valid_barcodes)

otu_melanie_long <- otu_melanie %>%
  rename_with(~ str_remove(.x, "e17a8f2887894f8d7becdbeaafbc97db14bc8e66_EXP-PBC096_")) %>%
  pivot_longer(-tax_id, names_to = "barcode", values_to = "count_melanie") %>%
  filter(barcode %in% valid_barcodes)

message("✅ Merged long-format OTU counts from Marly & Melanie.")

#----------------------------------------------------------*
# Step 4.1c: Load EMU Taxonomy Tables
#----------------------------------------------------------*
message("🔹 Loading taxonomy tables...")

tax_marly_path   <- file.path(processed_data, "EMU_output", "marly_standard_filtering", "taxonomy_marly_standard.tsv")
tax_melanie_path <- file.path(processed_data, "EMU_output", "melanie_lenient_filtering", "taxonomy_melanie_lenient.tsv")

tax_marly   <- read_tsv(tax_marly_path, show_col_types = FALSE)
tax_melanie <- read_tsv(tax_melanie_path, show_col_types = FALSE)

#----------------------------------------------------------*
# Step 4.1d: Merge OTU counts with taxonomy
#----------------------------------------------------------*
message("🔹 Annotating OTU tables with taxonomy...")

otu_marly_annotated <- otu_marly_long %>%
  left_join(tax_marly, by = "tax_id")

otu_melanie_annotated <- otu_melanie_long %>%
  left_join(tax_melanie, by = "tax_id")

#----------------------------------------------------------*
# Step 4.1e: Merge annotated OTUs with rodent metadata
#----------------------------------------------------------*
message("🔹 Merging annotated OTUs with rodent metadata...")

otu_marly_full <- otu_marly_annotated %>%
  left_join(rodent_data, by = "barcode")

otu_melanie_full <- otu_melanie_annotated %>%
  left_join(rodent_data, by = "barcode")

#----------------------------------------------------------*
# Step 4.1f: Save Output
#----------------------------------------------------------*
message("💾 Saving annotated OTU + taxonomy + metadata tables...")

write_csv(otu_marly_full,
          file.path(processed_data, "EMU_output", "marly_standard_filtering", "otu_taxonomy_metadata_marly.csv"))

write_csv(otu_melanie_full,
          file.path(processed_data, "EMU_output", "melanie_lenient_filtering", "otu_taxonomy_metadata_melanie.csv"))

# 🧹 Clean up
rm(
  otu_marly, otu_melanie,
  otu_marly_long, otu_melanie_long,
  otu_marly_annotated, otu_melanie_annotated,
  tax_marly, tax_melanie,
  otu_marly_path, otu_melanie_path
)

message("🎉 Done! OTU tables successfully merged and saved.")

