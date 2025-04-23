# ************************************************************
# Script: 04b_beta_pca_clr_taxa.R
# Purpose: PCA ordination with CLR-transformed taxa
# Input:   processed/phyloseq_16s_decontaminated.rds
# Output:  Results/Figures/16S_Beta/pca_clr_taxa.jpeg
# ************************************************************

message("🔬 Starting CLR-transformed PCA for genus-level taxa...")

# Load object
ps <- readRDS(file.path(processed_data, "phyloseq_16s_decontaminated.rds"))

# Aggregate to genus level
ps_genus <- tax_glom(ps, taxrank = "genus")

# CLR transform using microbiome package (returns a phyloseq object)
ps_clr <- microbiome::transform(ps_genus, "clr")

# Ordinate with PCA
pca_ord <- ordinate(ps_clr, method = "RDA")  # use RDA, not PCA, for CLR-transformed data

# Plot
p_pca <- plot_ordination(ps_clr, pca_ord, color = "sample_or_control") +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "PCA (CLR-transformed Genus)", color = "Sample Type")

# Save
ggsave(file.path(results_dir, "Figures", "16S_Beta", "pca_clr_taxa.jpeg"),
       plot = p_pca, width = 6, height = 4, dpi = 300)

######################### Cleaning ###########################################
# Clean temporary objects from your environment
rm(p_pca, pca_ord, ps, ps_clr, ps_genus)
