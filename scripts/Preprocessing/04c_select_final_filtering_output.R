# ***********************************************************
# Script: 04c_select_final_filtering_output.R
# Purpose: Finalize EMU Marly output and save for phyloseq input
# ***********************************************************

message("📌 Finalizing Marly's EMU 16S output for phyloseq...")

# -----------------------------------------------------------
# Step 1: Discard Melanie's OTU table
# -----------------------------------------------------------
rm(otu_melanie_full)

# -----------------------------------------------------------
# Step 2: Strip annotated table to essential OTU + taxonomy columns
# -----------------------------------------------------------
otu_stripped <- otu_marly_full %>%
  dplyr::select(
    barcode, tax_id, count_marly,
    species, genus, family, order, class, phylum, superkingdom
  )

# -----------------------------------------------------------
# Step 3: Merge with rodent metadata
# -----------------------------------------------------------
otu_joined <- otu_stripped %>%
  left_join(rodent_data, by = "barcode") %>%
  mutate(Gene = "16s")

# -----------------------------------------------------------
# Step 4: Rename for clarity
# -----------------------------------------------------------
rodent_data <- otu_joined %>%
  rename(
    barcode_16s = barcode,
    count_16s = count_marly
  )

# -----------------------------------------------------------
# Step 5: Create unique ID for sample-barcode combo
# -----------------------------------------------------------
rodent_data <- rodent_data %>%
  mutate(Sample_ID_phyloseq = paste(Sample_ID, barcode_16s, sep = "__"))

message("✅ Final merged dataset with 16S annotations and barcode IDs: ",
        nrow(rodent_data), " rows.")

# -----------------------------------------------------------
# Step 6: Export table for phyloseq
# -----------------------------------------------------------
otu_16s_phyloseq <- rodent_data %>%
  filter(!is.na(count_16s)) %>%
  dplyr::select(
    Sample_ID_phyloseq,  # <- unique sample-barcode combo
    tax_id,
    count_16s,
    species, genus, family, order, class, phylum, superkingdom,
    all_of(trapping_vars),
    conc_16s__PCR, Gene,
    Sample_ID, barcode_16s, is_control, sample_or_control  # keep traceability
  )

write_csv(
  otu_16s_phyloseq,
  file.path(processed_data, "EMU_output", "marly_standard_filtering", "otu_16s_for_phyloseq.csv")
)

message("💾 Saved cleaned 16S dataset for phyloseq construction.")

# -----------------------------------------------------------
# Step 7: Final cleanup – remove unneeded objects
# -----------------------------------------------------------
rm(
  barcode_ref,
  otu_stripped,
  otu_joined,
  otu_marly_full,
  otu_marly_path,
  otu_16s_phyloseq,
  tax_marly_path,
  tax_melanie_path,
  tax_marly,
  tax_melanie
)

# You may choose to keep `rodent_data`, `trapping_vars` if used downstream
message("🧹 Environment cleaned. Ready for phyloseq!")
