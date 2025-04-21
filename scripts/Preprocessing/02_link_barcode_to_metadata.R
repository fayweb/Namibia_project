# ***********************************************************
# Script: 02_link_barcode_to_metadata.R
# Title: Link Barcode IDs to Rodent Metadata
# Purpose: Merge barcode-to-sample reference with cleaned field metadata
# Author: Fay Webster
# ***********************************************************

message("🔗 Linking barcode reference data to cleaned rodent metadata...")

# Load barcode-sample reference file
barcode_ref <- read_csv("data/raw/16S_sequencing/reference_files/barcode_sample_reference.csv",
                        show_col_types = FALSE)

# 🧹 Clean barcode column: remove unnecessary prefix from ONT barcodes
barcode_ref <- barcode_ref %>%
  mutate(barcode = str_remove(barcode, "e17a8f2887894f8d7becdbeaafbc97db14bc8e66_EXP-PBC096_"))

# 🧼 Select only the clean columns you need and fix column names
barcode_ref <- barcode_ref %>%
  rename(Sample_ID = `Sample_ID...5`) %>%
  select(Sample_ID, barcode, conc_16s__PCR, sample_or_control)

# ✅ Preserve all sample-barcode pairings (even duplicates, blanks, etc.)
# Now join this to rodent metadata on Sample_ID (many-to-one allowed)

rodent_data <- full_join(barcode_ref, rodent_data, by = "Sample_ID")

# ✅ Flag any missing metadata (e.g. for blanks/controls) with a note
rodent_data <- rodent_data %>%
  mutate(is_control = if_else(sample_or_control == "control", TRUE, FALSE))

message("✅ Barcode metadata successfully linked to rodent data (rows: ", nrow(rodent_data), ")")

