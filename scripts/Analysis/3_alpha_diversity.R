# ************************************************************
# Script: 3_alpha_diversity.R
# Purpose: Generate rarefaction curves, summary plots, and alpha diversity comparisons
# Input:   processed/phyloseq_16s_decontaminated.rds
# Output:  Results/QC/16S/ (rarecurve, read depth histograms, diversity plots)
# ************************************************************

message("🌱 Starting alpha diversity analysis...")

# ------------------------------------------------------------
# Step 1: Load phyloseq object and apply basic filtering
# ------------------------------------------------------------
ps <- readRDS(file.path(processed_data, "phyloseq_16s_decontaminated.rds"))
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

read_summary <- summary(sample_sums(ps))
message("🧬 Read depth range: ", paste0(range(sample_sums(ps)), collapse = " - "))

# ------------------------------------------------------------
# Step 2: Generate Rarefaction Curve (Sample-wise)
# ------------------------------------------------------------
plot_rarefaction_with_filter <- function(physeq_obj, min_sum = 10, rarefy_to = 20000) {
  otu_tab <- t(as(otu_table(physeq_obj), "matrix"))
  metadata <- data.frame(sample_data(physeq_obj))

  otu_tab <- round(otu_tab)
  otu_tab[is.na(otu_tab)] <- 0
  otu_tab <- otu_tab[rowSums(otu_tab) >= min_sum, ]
  metadata <- metadata[rownames(otu_tab), ]

  if (nrow(otu_tab) < 2) {
    warning("Not enough samples after filtering for rarefaction.")
    return(NULL)
  }

  rare_depth <- min(min(rowSums(otu_tab)), rarefy_to)
  group_colors <- if ("sample_or_control" %in% colnames(metadata)) {
    cols <- c("sample" = "black", "control" = "red")
    cols[metadata$sample_or_control]
  } else {
    rainbow(nrow(otu_tab))
  }

  jpeg(file.path(results_dir, "Figures", "16S_QC", "rarefaction_curve.jpeg"), width = 8, height = 5, units = "in", res = 300)
  vegan::rarecurve(otu_tab, step = 100, label = FALSE, sample = rare_depth, col = group_colors, cex = 0.6, main = "Rarefaction Curve by Sample Type")
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
# Step 4: Plot Median DNA Concentration
# ------------------------------------------------------------
p_dna_conc <- ggplot(sample_dat_info, aes(x = sample_or_control, y = conc_16s__PCR)) +
  geom_boxplot(aes(fill = sample_or_control), alpha = 0.6, outlier.shape = NA, width = 0.5) +
  theme_minimal() +
  labs(x = "Sample Type", y = "DNA Concentration (PCR)", title = "DNA Concentration by Sample Type") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14), legend.position = "none")

ggsave(filename = file.path(results_dir, "Figures", "16S_QC", "dna_concentration_by_sample.jpeg"),
       plot = p_dna_conc, width = 6, height = 4, dpi = 300)

# ************************************************************
# Alpha diversity (Shannon, Observed, InvSimpson)
# ************************************************************
message("📈 Calculating alpha diversity (Shannon, Observed, InvSimpson)...")

otu_table(ps) <- otu_table(round(otu_table(ps)), taxa_are_rows = TRUE)

alpha_df <- phyloseq::estimate_richness(ps, measures = c("Observed", "Shannon", "InvSimpson")) %>%
  tibble::rownames_to_column("Sample_ID") %>%
  left_join(sample_data(ps) %>% data.frame() %>% tibble::rownames_to_column("Sample_ID"), by = "Sample_ID")

# Create folder for saving alpha diversity table
dir.create(file.path(results_dir, "Tables", "16S_QC"), recursive = TRUE, showWarnings = FALSE)

# Save alpha diversity metrics
write_csv(alpha_df, file.path(results_dir, "Tables", "16S_QC", "alpha_diversity_metrics.csv"))

# Panel plot: diversity metrics per sample type
p_alpha <- alpha_df %>%
  pivot_longer(cols = c(Observed, Shannon, InvSimpson), names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = sample_or_control, y = Value, fill = sample_or_control)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
  facet_wrap(~Metric, scales = "free_y") +
  theme_minimal() +
  labs(title = "Alpha Diversity by Sample Type", x = "Sample Type", y = "Diversity Metric Value") +
  theme(strip.text = element_text(size = 13, face = "bold"),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        legend.position = "none")

ggsave(file.path(results_dir, "Figures", "16S_QC", "alpha_diversity_panel.jpeg"),
       plot = p_alpha, width = 8, height = 5, dpi = 300)

message("📈 Saved alpha diversity panel to: Figures/16S_QC/alpha_diversity_panel.jpeg")

# ************************************************************
# Compare Alpha Diversity by Rodent Species
# ************************************************************

message("📈 Creating alpha diversity boxplot (ordered by Shannon diversity)...")

# Reorder species by median Shannon diversity
alpha_df_ordered <- alpha_df %>%
  group_by(Morphology_species) %>%
  summarise(median_Shannon = median(Shannon, na.rm = TRUE)) %>%
  arrange(median_Shannon)

alpha_df <- alpha_df %>%
  mutate(Morphology_species = factor(Morphology_species, levels = alpha_df_ordered$Morphology_species))

# Boxplot for alpha diversity by species
p_alpha_species <- ggplot(alpha_df, aes(x = Morphology_species, y = Shannon, fill = Morphology_species)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  theme_minimal() +
  labs(title = "Alpha Diversity (Shannon) by Species", x = "Rodent Species", y = "Shannon Diversity Index") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

ggsave(file.path(results_dir, "Figures", "16S_QC", "alpha_diversity_species_ordered.jpeg"),
       plot = p_alpha_species, width = 8, height = 5, dpi = 300)

message("✅ Alpha diversity by species (ordered) saved.")

# ************************************************************
# Compare Alpha Diversity by Trapping Location
# ************************************************************

message("📈 Creating alpha diversity boxplot (Shannon) by Trapping Location...")

# Boxplot for Shannon diversity by trapping location
p_alpha_location <- ggplot(alpha_df, aes(x = Location_type, y = Shannon, fill = Location_type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
  theme_minimal(base_size = 15) +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Alpha Diversity (Shannon) by Trapping Location", x = "Trapping Location", y = "Shannon Diversity Index") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14), legend.position = "none", plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

ggsave(file.path(results_dir, "Figures", "16S_QC", "alpha_diversity_location.jpeg"),
       plot = p_alpha_location, width = 8, height = 5, dpi = 300)

message("✅ Alpha diversity by location (Shannon) saved.")

# ************************************************************
# Final Step: Clean up temporary objects
# ************************************************************

rm(ps, read_summary, sample_dat_info, dna_conc_summary, alpha_df_ordered,
   p_dna_conc, p_alpha_species, p_alpha_location, alpha_df, p_alpha)

message("🧹 Cleaned up temporary objects.")

