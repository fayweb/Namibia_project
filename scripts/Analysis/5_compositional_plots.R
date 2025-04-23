# ************************************************************
# Script: 5_compositional_plots.R
# Purpose: Create compositional plots for relative abundance
# Input:   processed/phyloseq_16s_decontaminated.rds
# Output:  Results/Figures/16S_QC/compositional_plots.jpeg
# ************************************************************

message("🔬 Starting the creation of compositional plots...")

# ------------------------------------------------------------
# Step 1: Load phyloseq object and apply filters
# ------------------------------------------------------------
ps <- readRDS(file.path(processed_data, "phyloseq_16s_decontaminated.rds"))

# Filter samples with fewer than 200 reads
ps <- prune_samples(sample_sums(ps) >= 200, ps)

# ------------------------------------------------------------
# Step 2: Clean the taxonomy table
# ------------------------------------------------------------
# Convert tax_table to a data frame for easier manipulation
tax_df <- as.data.frame(tax_table(ps))

# Remove rows with missing species or genus
tax_df <- tax_df[!is.na(tax_df$species) & !is.na(tax_df$genus), ]

# Replace missing species/genus with "Unknown"
tax_df$species[is.na(tax_df$species)] <- "Unknown"
tax_df$genus[is.na(tax_df$genus)] <- "Unknown"

# Remove unwanted characters like semicolons and UCG terms
tax_df$species <- gsub("UCG.*", "", tax_df$species)  # Remove UCG terms
tax_df$species <- gsub(";", "", tax_df$species)  # Remove semicolons

# Convert back to matrix and update tax_table
tax_table(ps) <- as.matrix(tax_df)

# ------------------------------------------------------------
# Step 3: Apply tax_fix() to standardize the taxonomy
# ------------------------------------------------------------
ps2 <- ps %>%
  tax_fix(min_length = 4, unknowns = c(""), sep = " ", anon_unique = TRUE, suffix_rank = "classified") %>%
  phyloseq_validate()

# ------------------------------------------------------------
# Step 4: Handle problematic genera using tax_fix
# ------------------------------------------------------------
# Apply tax_fix to resolve convergent taxa (e.g., Lactobacillus, Bacteroides)
ps2 <- ps2 %>%
  tax_fix(unknowns = c("Alistipes", "Aquitalea", "Asaia", "Bacillus", "Bacteroides",
                       "Bifidobacterium", "Deinococcus", "Enterococcus", "Geodermatophilus",
                       "Lactobacillus", "Neisseria", "Parabacteroides", "Paracoccus", "Salmonella")) %>%
  phyloseq_validate()

cleaned_ps2_path <- file.path(processed_data, "phyloseq_16s_cleaned_genera.rds")
saveRDS(ps2, file = cleaned_ps2_path)

# ------------------------------------------------------------
# Step 5: Create compositional barplot at the Genus level
# ------------------------------------------------------------
if ("genus" %in% rank_names(ps2)) {
  p_compositional <- ps2 %>%
    comp_barplot(tax_level = "genus") +
    coord_flip()  # Horizontal bars for better readability

  # Save compositional plot
  ggsave(file.path(results_dir, "Figures", "16S_QC", "compositional_plots.jpeg"),
         plot = p_compositional, width = 8, height = 5)

  message("✅ Compositional plot saved.")
} else {
  message("🚨 'Genus' rank is not available in the taxonomy table. Please check the dataset.")
}

# ------------------------------------------------------------

# ------------------------------------------------------------
# Final Step: Clean up temporary objects
# ------------------------------------------------------------

# Remove unnecessary objects from the environment to free up memory
rm(list = c("ps", "tax_df", "p_compositional", "ps2"))

# Optional: If there are any other temporary objects you used that are not needed further
# You can add them to the list like this:
# rm(list = c("other_temp_object_1", "other_temp_object_2"))

message("🧹 Cleaned up temporary objects.")

