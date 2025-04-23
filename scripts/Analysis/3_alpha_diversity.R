# ************************************************************
# Script: 3_alpha_diversity.R
# Purpose: Generate rarefaction curves and summary plots for 16S alpha diversity
# Input:   processed/phyloseq_16s_decontaminated.rds
# Output:  Results/QC/16S/ (rarecurve, read depth histograms)
# ************************************************************

message("🌱 Starting alpha diversity analysis...")

# ------------------------------------------------------------
# Step 1: Load phyloseq object and basic filtering
# ------------------------------------------------------------
ps <- readRDS(file.path(processed_data, "phyloseq_16s_decontaminated.rds"))
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

read_summary <- summary(sample_sums(ps))
message("🧬 Read depth range: ", paste0(range(sample_sums(ps)),

                                        collapse = " - "))

# ------------------------------------------------------------
# Step 2: Rarefaction Curve (Sample-wise)
# ------------------------------------------------------------
plot_rarefaction_with_filter <- function(physeq_obj, min_sum = 10,
                                         rarefy_to = 20000) {
  otu_tab <- t(as(otu_table(physeq_obj), "matrix"))
  metadata <- data.frame(sample_data(physeq_obj))

  # Round to integer counts
  otu_tab <- round(otu_tab)

  # Replace NAs with 0 and filter
  otu_tab[is.na(otu_tab)] <- 0
  otu_tab <- otu_tab[rowSums(otu_tab) >= min_sum, ]
  metadata <- metadata[rownames(otu_tab), ]

  if (nrow(otu_tab) < 2) {
    warning("Not enough samples after filtering for rarefaction.")
    return(NULL)
  }

  # Cap rarefaction depth to a safe minimum
  rare_depth <- min(min(rowSums(otu_tab)), rarefy_to)

  # Assign color
  group_colors <- if ("sample_or_control" %in% colnames(metadata)) {
    cols <- c("sample" = "black", "control" = "red")
    cols[metadata$sample_or_control]
  } else {
    rainbow(nrow(otu_tab))
  }

  # Save plot
  jpeg(file.path(results_dir, "Figures", "16S_QC", "rarefaction_curve.jpeg"),
       width = 8, height = 5, units = "in", res = 300)

  vegan::rarecurve(otu_tab, step = 100, label = FALSE,
                   sample = rare_depth, col = group_colors,
                   cex = 0.6, main = "Rarefaction Curve by Sample Type")

  dev.off()
  message("📈 Rarefaction curve saved successfully.")
}

plot_rarefaction_with_filter(ps)



# ------------------------------------------------------------
# Step 3: DNA Concentration Summary
# ------------------------------------------------------------
sample_dat_info <- data.frame(sample_data(ps))

dna_conc_summary <- sample_dat_info %>%
  group_by(sample_or_control) %>%
  summarise(
    median_conc = median(conc_16s__PCR, na.rm = TRUE),
    sd_conc = sd(conc_16s__PCR, na.rm = TRUE),
    n = n()
  )

print(dna_conc_summary)

# ------------------------------------------------------------
# Step 4: Plot median DNA concentration
# ------------------------------------------------------------
p_dna_conc <- ggplot(sample_dat_info, aes(x = sample_or_control,
                                          y = conc_16s__PCR)) +
  geom_boxplot(aes(fill = sample_or_control), alpha = 0.6,
               outlier.shape = NA, width = 0.5) +
  theme_minimal() +
  labs(x = "Sample Type", y = "DNA Concentration (PCR)",
       title = "DNA Concentration by Sample Type") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "none")

ggsave(filename = file.path(results_dir, "Figures", "16S_QC",
                            "dna_concentration_by_sample.jpeg"),
       plot = p_dna_conc, width = 6, height = 4, dpi = 300)

# ************************************************************
# Alpha diversity
# ************************************************************

message("📈 Calculating alpha diversity (Shannon, Observed, InvSimpson)...")

# ------------------------------------------------------------
# Step 3: Visualize Alpha Diversity (Shannon + Others)
# ------------------------------------------------------------

# Round OTU table to integers (required by estimate_richness)
otu_table(ps) <- otu_table(round(otu_table(ps)), taxa_are_rows = TRUE)

# Compute alpha diversity metrics
alpha_df <- phyloseq::estimate_richness(ps,
                                        measures =
                                          c("Observed", "Shannon", "InvSimpson")) %>%
  tibble::rownames_to_column("Sample_ID") %>%
  left_join(sample_data(ps) %>%
              data.frame() %>%
              tibble::rownames_to_column("Sample_ID"),
            by = "Sample_ID")

# Note: estimate_richness warning about missing singletons is expected
# because the dataset has been decontaminated and filtered
# Create folder if it doesn't exist
dir.create(file.path(results_dir, "Tables", "16S_QC"),
           recursive = TRUE, showWarnings = FALSE)

# Save alpha diversity table
write_csv(alpha_df, file.path(results_dir, "Tables", "16S_QC",
                              "alpha_diversity_metrics.csv"))


# Panel plot: diversity metrics per sample type
p_alpha <- alpha_df %>%
  pivot_longer(cols = c(Observed, Shannon, InvSimpson),
               names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = sample_or_control, y = Value, fill = sample_or_control)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
  facet_wrap(~Metric, scales = "free_y") +
  theme_minimal() +
  labs(title = "Alpha Diversity by Sample Type",
       x = "Sample Type",
       y = "Diversity Metric Value") +
  theme(
    strip.text = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    legend.position = "none"
  )

# Save plot
ggsave(file.path(results_dir, "Figures", "16S_QC",
                 "alpha_diversity_panel.jpeg"),
       plot = p_alpha, width = 8, height = 5, dpi = 300)

message(
  "📈 Saved alpha diversity panel to: Figures/16S_QC/alpha_diversity_panel.jpeg")


# ------------------------------------------------------------
# Step 5: Visualize Alpha Diversity (Shannon + Others)
# ------------------------------------------------------------

# Round OTU table to integers (required by estimate_richness)
otu_table(ps) <- otu_table(round(otu_table(ps)), taxa_are_rows = TRUE)

# Compute alpha diversity metrics
alpha_df <- phyloseq::estimate_richness(
  ps, measures =
    c("Observed", "Shannon", "InvSimpson")) %>%
  tibble::rownames_to_column("Sample_ID") %>%
  left_join(sample_data(ps) %>%
              data.frame() %>%
              tibble::rownames_to_column("Sample_ID"),
            by = "Sample_ID")

# Note: estimate_richness warning about missing singletons is expected
# because the dataset has been decontaminated and filtered
# Create folder if it doesn't exist
dir.create(file.path(results_dir, "Tables", "16S_QC"),
           recursive = TRUE, showWarnings = FALSE)

# Save alpha diversity table
write_csv(alpha_df, file.path(results_dir, "Tables", "16S_QC",
                              "alpha_diversity_metrics.csv"))

# ************************************************************
# Compare alpha diversity by Rodent Species and Trapping Location
# ************************************************************

# Boxplot: Diversity metrics per species
p_alpha_species <- alpha_df %>%
  ggplot(aes(x = Morphology_species, y = Shannon, fill = Morphology_species)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  theme_minimal() +
  labs(title = "Alpha Diversity (Shannon) by Species",
       x = "Rodent Species", y = "Shannon Diversity Index") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save species diversity plot
ggsave(file.path(results_dir, "Figures", "16S_QC",
                 "alpha_diversity_species.jpeg"),
       plot = p_alpha_species, width = 8, height = 5, dpi = 300)

# Boxplot: Diversity metrics per trapping location
p_alpha_location <- alpha_df %>%
  ggplot(aes(x = Location_type, y = Shannon, fill = Location_type)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  theme_minimal() +
  labs(title = "Alpha Diversity (Shannon) by Trapping Location",
       x = "Trapping Location", y = "Shannon Diversity Index") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save location diversity plot
ggsave(file.path(results_dir, "Figures", "16S_QC",
                 "alpha_diversity_location.jpeg"),
       plot = p_alpha_location, width = 8, height = 5, dpi = 300)

message("📈 Saved alpha diversity comparison plots by species and location.")

# ------------------------------------------------------------
# Step 6: Clean environment
## ------------------------------------------------------------
rm(ps, read_summary, sample_dat_info, dna_conc_summary,
   p_dna_conc, p_alpha_species, p_alpha_location, alpha_df)

message("🧹 Cleaned up temporary objects.")
