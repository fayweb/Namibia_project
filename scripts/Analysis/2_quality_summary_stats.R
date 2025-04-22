# ************************************************************
# Script: 2a_quality_summary_stats.R
# Purpose: Summarize read depth, abundance, and prevalence distributions
# Input:   processed/phyloseq_16s_decontaminated.rds
# Output:  Results/QC/16S/ (various plots)
# ************************************************************

message("📊 Summarizing 16S quality metrics (library size, abundance, prevalence)...")

# ------------------------------------------------------------
# Step 1: Load phyloseq object
# ------------------------------------------------------------
ps <- readRDS(file.path(processed_data, "phyloseq_16s_decontaminated.rds"))

# ------------------------------------------------------------
# Step 2: Library Size Plot by Sample Type
# ------------------------------------------------------------
sample_data(ps)$LibrarySize <- sample_sums(ps)

p_libsize <- ggplot(as.data.frame(sample_data(ps)), aes(
  x = reorder(rownames(sample_data(ps)), LibrarySize),
  y = LibrarySize,
  color = sample_or_control)) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  labs(
    x = "Sample Index",
    y = "Library Size (reads)",
    title = "Read Depth by Sample Type"
  )

ggsave(file.path(results_dir, "Figures", "16S_QC", "library_sizes_by_type.jpeg"),
       plot = p_libsize, width = 8, height = 5)

# ------------------------------------------------------------
# Step 3: Histogram of Raw Abundances (log scale)
# ------------------------------------------------------------
ps_df_taxa <- data.table(
  tax_table(ps),
  ASVabundance = taxa_sums(ps),
  ASV = taxa_names(ps)
)

p_abund_hist <- ggplot(ps_df_taxa, aes(ASVabundance)) +
  geom_histogram() +
  scale_x_log10() +
  theme_bw() +
  labs(title = "Histogram of Raw Taxa Counts (log10)", y = "Frequency", x = "ASV Abundance")

ggsave(file.path(results_dir,  "Figures", "16S_QC", "histogram_raw_abundance.jpeg"),
       plot = p_abund_hist, width = 6, height = 4)

# ------------------------------------------------------------
# Step 4: Coefficient of Variation vs. Mean Abundance
# ------------------------------------------------------------

# 1. Convert to relative abundance
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# 2. Extract OTU table
otu_rel_df <- as.data.frame(otu_table(ps_rel))

# 3. Transpose so taxa are rows
if (!taxa_are_rows(ps_rel)) {
  otu_rel_df <- t(otu_rel_df)
}

# 4. Compute mean and sd per taxon
cv_df <- as.data.frame(otu_rel_df) %>%
  rownames_to_column("taxon") %>%
  mutate(
    mean_abundance = rowMeans(select(., -taxon), na.rm = TRUE),
    sd_abundance = apply(select(., -taxon), 1, sd, na.rm = TRUE),
    cv = sd_abundance / mean_abundance
  )

# 5. Plot: CV vs mean abundance
p_cv <- ggplot(cv_df, aes(x = mean_abundance, y = cv)) +
  geom_point(alpha = 0.6) +
  scale_x_log10() +
  scale_y_log10() +
  theme_minimal() +
  labs(
    title = "Coefficient of Variation vs Mean Abundance",
    x = "Mean Relative Abundance (log10)",
    y = "Coefficient of Variation (log10)",
  ) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

print(p_cv)

# Save
ggsave(filename = file.path(results_dir, "Figures", "16S_QC", "cv_vs_mean_abundance.jpeg"),
       plot = p_cv, width = 6, height = 4, dpi = 300)

# ------------------------------------------------------------
# Step 5: Prevalence by Phylum
# ------------------------------------------------------------

# Step 1: Prepare prevalence dataframe manually
tax_df <- as.data.frame(tax_table(ps))
otu_mat <- as.matrix(otu_table(ps))

if (!taxa_are_rows(ps)) {
  otu_mat <- t(otu_mat)
}

# Calculate prevalence for each taxon (number of samples where present)
tax_df$prevalence <- rowSums(otu_mat > 0)
tax_df$phylum <- tax_df$phylum %>% replace_na("Unclassified")

# Step 2: Aggregate prevalence by phylum
prev_summary <- tax_df %>%
  group_by(phylum) %>%
  summarise(
    prevalence_mean = mean(prevalence),
    prevalence_sd = sd(prevalence),
    taxa_count = n()
  ) %>%
  arrange(desc(prevalence_mean))

# Step 3: Plot
p_prev <- ggplot(prev_summary, aes(x = reorder(phylum, -prevalence_mean),
                                   y = prevalence_mean)) +
  geom_col(fill = "deeppink") +
  geom_errorbar(aes(ymin = prevalence_mean - prevalence_sd,
                    ymax = prevalence_mean + prevalence_sd),
                width = 0.2, color = "grey40") +
  theme_minimal() +
  labs(
    title = "Average Taxon Prevalence by Phylum",
    x = "Phylum",
    y = "Average Prevalence Across Samples"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

# Save plot
ggsave(filename = file.path(results_dir, "Figures", "16S_QC", "prevalence_by_phylum.jpeg"),
       plot = p_prev, width = 8, height = 5, dpi = 300)


# ------------------------------------------------------------
# Step 6: Rare Taxa Summary
# ------------------------------------------------------------

# Total number of taxa
ntaxa_total <- ntaxa(ps)

# Subset taxa with known phylum
ps_present <- subset_taxa(ps, !is.na(tax_table(ps)[, "phylum"])
                          & tax_table(ps)[, "phylum"] != "")

# Total reads from those taxa
reads_total <- sum(taxa_sums(ps_present))

# How many taxa have fewer than 10 reads
low_abundance_taxa <- taxa_sums(ps_present) < 10
n_low_abund <- sum(low_abundance_taxa)
reads_low_abund <- sum(taxa_sums(ps_present)[low_abundance_taxa])

# What proportion is this?
prop_low_abund <- round((n_low_abund / ntaxa(ps_present)) * 100, 2)

# Singleton and doubleton counts
n_singletons <- sum(taxa_sums(ps_present) == 1)
n_doubletons <- sum(taxa_sums(ps_present) == 2)
n_single_double <- n_singletons + n_doubletons
prop_single_double <- round((n_single_double / ntaxa(ps_present)) * 100, 2)

# Print summary
message("📊 Total taxa (with phylum): ", ntaxa(ps_present))
message("📉 Taxa with <10 reads: ", n_low_abund, " (", prop_low_abund, "%)")
message("📦 Total reads: ", reads_total)
message("🧬 Reads in rare taxa (<10): ", reads_low_abund)
message("🧫 Singletons: ", n_singletons, " | Doubletons: ", n_doubletons)
message("⚠️ Singleton + doubleton = ", prop_single_double, "% of phylum-assigned taxa")

# ------------------------------------------------------------
# Step 7: Clean up
# ------------------------------------------------------------
rm(
  ps, ps_rel, ps_phylum, ps_present,
  tax_df, otu_mat, otu_rel_df, cv_df, prev_summary,
  low_abundance_taxa,
  p_libsize, p_abund_hist, p_cv, p_prev
)

message("🧹 Cleaned up temporary objects.")


# ------------------------------------------------------------
# Step 7: Clean up
# ------------------------------------------------------------
rm(ps, p_libsize, p_abund_hist, p_cv, p_prev,
   ps_df_taxa, ps_present,
   rare_counts, rare_singletons, rare_doubletons,
   rare_lessthan10, reads_lessthan10)

message("🧹 Cleaned up temporary objects.")
