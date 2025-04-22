# ************************************************************
# Script: 04a_beta_distance_ordination.R
# Purpose: Beta diversity using Bray-Curtis & NMDS
# Input:   processed/phyloseq_16s_decontaminated.rds
# Output:  Results/Figures/16S_Beta/nmds_sample_type.jpeg
# ************************************************************

message("🧭 Starting beta diversity ordination...")

# Load object
ps <- readRDS(file.path(processed_data, "phyloseq_16s_decontaminated.rds"))

# Transform to relative abundance
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# Compute Bray-Curtis distance
bray_dist <- phyloseq::distance(ps_rel, method = "bray")

# NMDS ordination
nmds_ord <- ordinate(ps_rel, method = "NMDS", distance = bray_dist)

# Plot NMDS
p_nmds <- plot_ordination(ps_rel, nmds_ord, color = "sample_or_control") +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "NMDS Ordination by Sample Type", color = "Sample Type")

# Save
ggsave(file.path(results_dir, "Figures", "16S_Beta", "nmds_sample_type.jpeg"),
       plot = p_nmds, width = 6, height = 4, dpi = 300)

message("✅ NMDS ordination saved.")
rm(ps, ps_rel, bray_dist, nmds_ord, p_nmds)
