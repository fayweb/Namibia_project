# 04c_select_final_filtering_output.R
# Purpose: Merge selected EMU output (Marly) with full rodent metadata
# we are removing Melanie's lenient filtering strategy data
rm(otu_melanie_full)

# 1. First keep original name so we can join
otu_stripped <- otu_marly_full %>%
  select(barcode, tax_id, count_marly, species, genus, family, order, class,
         phylum, superkingdom)

# 2. Join with rodent_data first (barcode is still available for matching)
otu_joined <- otu_stripped %>%
  left_join(rodent_data, by = "barcode")

Gene <- "16s"

otu_joined <- cbind(otu_joined, Gene)

# 3. Then rename columns for clarity
rodent_data <- otu_joined %>%
  rename(
    barcode_16s = barcode,
    count_16s = count_marly
    # Keep tax_id unsuffixed for unified use
  )


# 5. Sanity check
message("✅ Final merged dataset with 16S columns ready: ",
        nrow(rodent_data), " rows.")


rm(otu_joined, otu_stripped, otu_marly_full, barcode_ref)
