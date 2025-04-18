# ***********************************************************
# Script: 01_construct_phyloseq.R
# Purpose: Construct a phyloseq object from preprocessed 16S data
# ***********************************************************
# Script: 01_construct_phyloseq.R
# Purpose: Construct a phyloseq object from Marly's 16S OTU dataset
# Input:   data/processed/EMU_output/marly_standard_filtering/otu_16s_for_phyloseq.csv
# Output:  phyloseq_16s_final.rds (to be used for alpha/beta/diversity plots)
# ***********************************************************

message("📦 Constructing phyloseq object from preprocessed 16S dataset...")

# -----------------------------------------------------------
# Step 1: Load preprocessed OTU + taxonomy + metadata table
# -----------------------------------------------------------
otu_path <- file.path(processed_data, "EMU_output", "marly_standard_filtering", "otu_16s_for_phyloseq.csv")
otu_df <- readr::read_csv(otu_path, show_col_types = FALSE)

message("✔️ Loaded: ", nrow(otu_df), " rows of 16S data with metadata")

# -----------------------------------------------------------
# Step 2: Create OTU table (rows = taxa, cols = samples)
# -----------------------------------------------------------
otu_long <- otu_df %>%
  dplyr::select(Sample_ID, tax_id, count_16s)

# -----------------------------------------------------------
# Step 2: Create OTU table (rows = taxa, cols = samples)
# -----------------------------------------------------------
otu_wide <- otu_long %>%
  dplyr::group_by(tax_id, Sample_ID) %>%
  dplyr::summarise(count_16s = sum(count_16s, na.rm = TRUE),
                   .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = Sample_ID,
    values_from = count_16s,
    values_fill = list(count_16s = 0)
  ) %>%
  tibble::column_to_rownames("tax_id")

otu_table_ps <- phyloseq::otu_table(as.matrix(otu_wide), taxa_are_rows = TRUE)


# -----------------------------------------------------------
# Step 3: Create Taxonomy Table (rows = taxa)
# -----------------------------------------------------------
taxonomy_df <- otu_df %>%
  dplyr::select(tax_id, species, genus, family, order, class, phylum, superkingdom) %>%
  dplyr::distinct(tax_id, .keep_all = TRUE) %>%
  tibble::column_to_rownames("tax_id")

tax_table_ps <- phyloseq::tax_table(as.matrix(taxonomy_df))

# -----------------------------------------------------------
# Step 4: Create Sample Metadata Table (rows = samples)
# -----------------------------------------------------------
sample_metadata_df <- otu_df %>%
  dplyr::select(Sample_ID, all_of(trapping_vars), conc_16s__PCR, Gene) %>%
  dplyr::distinct(Sample_ID, .keep_all = TRUE) %>%
  dplyr::filter(!is.na(Sample_ID)) %>%
  tibble::column_to_rownames("Sample_ID")

sample_data_ps <- phyloseq::sample_data(sample_metadata_df)

# -----------------------------------------------------------
# Step 5: Construct phyloseq object
# -----------------------------------------------------------
ps_16s_final <- phyloseq::phyloseq(
  otu_table_ps,
  tax_table_ps,
  sample_data_ps
)

message("✅ Successfully constructed phyloseq object with ",
        phyloseq::nsamples(ps_16s_final), " samples and ",
        phyloseq::ntaxa(ps_16s_final), " taxa.")

# -----------------------------------------------------------
# Step 6: Save output
# -----------------------------------------------------------
phyloseq_out_path <- file.path(processed_data, "phyloseq_16s_final.rds")
saveRDS(ps_16s_final, file = phyloseq_out_path)

message("💾 Saved phyloseq object to: ", phyloseq_out_path)

# -----------------------------------------------------------
# Step 7: Clean Environment (Optional but Recommended)
# -----------------------------------------------------------
# Remove intermediate objects used only for construction
rm(
  otu_df, otu_long, otu_wide,
  taxonomy_df, tax_table_ps,
  sample_metadata_df, sample_data_ps
)

message("🧹 Cleaned up temporary objects from environment.")

