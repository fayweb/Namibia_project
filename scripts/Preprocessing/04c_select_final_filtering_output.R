# 04c_select_final_filtering_output.R
# Purpose: Merge selected EMU output (Marly) with full rodent metadata
# we are removing Melanie's lenient filtering strategy data
rm(otu_melanie_full)

# 1. Strip to counts + taxonomy
otu_stripped <- otu_marly_full %>%
  select(barcode, tax_id, count_marly, species, genus, family, order, class,
         phylum, superkingdom)

otu_renamed <- otu_stripped %>%
  rename_with(.fn = ~ paste0(.x, "_16s"), .cols = -c(barcode, tax_id))


# 3. Join with rodent metadata
rodent_data <- otu_renamed %>%
  left_join(rodent_data, by = "barcode")

# 5. Sanity check
message("✅ Final merged dataset with 16S columns ready: ",
        nrow(rodent_data), " rows.")


rm(otu_renamed, otu_stripped, otu_marly_full)
