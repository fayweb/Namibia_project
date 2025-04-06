# ***********************************************************
# Title: Link Barcode IDs to Rodent Metadata
# Purpose: Merge barcode-to-sample reference with cleaned field metadata
# ***********************************************************
# Load barcode-sample reference file
barcode_ref <- read_csv("data/raw/16S_sequencing/reference_files/barcode_sample_reference.csv")

# remove the unecessary string in front of the barcode names
barcode_ref$barcode <- str_remove(
  string = barcode_ref$barcode,
  pattern = "e17a8f2887894f8d7becdbeaafbc97db14bc8e66_EXP-PBC096_")

# remove repate sample name column
barcode_ref <- barcode_ref %>%
  dplyr::select(-Sample_ID...6)

# rename now barcode_ref column sample id
barcode_ref <- barcode_ref %>%
  dplyr::rename(Sample_ID = Sample_ID...5) %>%
  dplyr::select(Sample_ID, conc_16s__PCR, barcode)

# Merge by Sample_ID
rodent_data <- left_join(rodent_data, barcode_ref, by = "Sample_ID")

message("✅ Barcode metadata linked to rodent data ")
