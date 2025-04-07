# ***********************************************************
# Title: Merge OTU Counts from EMU to Rodent Metadata
# Purpose: Integrate Marly's and Melanie's filtered OTU tables into rodent_data
# Author: Fay Webster
# ***********************************************************
# Purpose: Join OTU counts from Marly's and Melanie's EMU outputs
# into a single long-format table for easy comparison.
# This includes:
#   - tax_id (unique taxa)
#   - barcode (sequencing ID)
#   - count_marly (standard filtering)
#   - count_melanie (lenient filtering)
# ----------------------------------------------------------*

# Paths to OTU files
otu_marly_path <- "data/processed/EMU_output/marly_standard_filtering/otu_counts_marly_standard.tsv"
otu_melanie_path <- "data/processed/EMU_output/melanie_lenient_filtering/otu_counts_melanie_lenient.tsv"

# Load and pivot to long format
otu_marly <- read_tsv(otu_marly_path)
otu_melanie <- read_tsv(otu_melanie_path)

# Load barcode metadata (already cleaned and integrated into rodent_data)
barcode_vector <- unique(rodent_data$barcode)

# Clean and pivot Marly's OTU table
otu_marly_long <- otu_marly %>%
  rename_with(~ str_remove(.x, "e17a8f2887894f8d7becdbeaafbc97db14bc8e66_EXP-PBC096_")) %>%
  pivot_longer(-tax_id, names_to = "barcode", values_to = "count_marly") %>%
  filter(barcode %in% barcode_vector)

# Clean and pivot Melanie's OTU table
otu_melanie_long <- otu_melanie %>%
  rename_with(~ str_remove(.x, "e17a8f2887894f8d7becdbeaafbc97db14bc8e66_EXP-PBC096_")) %>%
  pivot_longer(-tax_id, names_to = "barcode", values_to = "count_melanie") %>%
  filter(barcode %in% barcode_vector)

otu_merged <- full_join(otu_marly_long, otu_melanie_long, by = c("tax_id", "barcode"))

# remove unecessary files
rm(otu_marly_path, otu_melanie_path, otu_marly, otu_melanie, otu_marly_long,
   otu_melanie_long)

# Save for reference (optional)
write_csv(otu_merged, file.path(processed_data, "EMU_output", "otu_counts_long_combined.csv"))
message("✅ Merged long-format OTU counts from Marly & Melanie saved to: data/processed/EMU_output/otu_counts_long_combined.csv")
