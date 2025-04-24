# ************************************************************
# Script: 04a_beta_distance_ordination.R
# Purpose: Beta diversity using Bray-Curtis & NMDS
# Input:   processed/phyloseq_16s_decontaminated.rds
# Output:  Results/Figures/16S_Beta/nmds_sample_type.jpeg
# ************************************************************

message("🧭 Starting beta diversity ordination...")

# ------------------------------------------------------------
# Step 1: Load the phyloseq object and transform to relative abundance
# ------------------------------------------------------------
ps <- readRDS(file.path(processed_data, "phyloseq_16s_decontaminated.rds"))

# Transform OTU table to relative abundance
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# ------------------------------------------------------------
# Step 2: Compute Bray-Curtis distance
# ------------------------------------------------------------
bray_dist <- phyloseq::distance(ps_rel, method = "bray")

# ------------------------------------------------------------
# Step 3: Perform NMDS ordination
# ------------------------------------------------------------
nmds_ord <- ordinate(ps_rel, method = "NMDS", distance = bray_dist)

# ------------------------------------------------------------
# Step 4: Plot NMDS ordination by sample type (default)
# ------------------------------------------------------------
p_nmds_sample_type <- plot_ordination(ps_rel, nmds_ord, color = "sample_or_control") +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "NMDS Ordination by Sample Type", color = "Sample Type")

# Save NMDS plot by sample type
ggsave(file.path(results_dir, "Figures", "16S_Beta", "nmds_sample_type.jpeg"),
       plot = p_nmds_sample_type, width = 6, height = 4, dpi = 300)

message("✅ NMDS ordination by sample type saved.")

# ------------------------------------------------------------
# Step 5: Explore Beta Diversity by Morphology Species
# ------------------------------------------------------------
p_nmds_species <- plot_ordination(ps_rel, nmds_ord, color = "Morphology_species") +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "NMDS Ordination by Rodent Species", color = "Rodent Species")

# Save NMDS plot by species
ggsave(file.path(results_dir, "Figures", "16S_Beta", "nmds_species.jpeg"),
       plot = p_nmds_species, width = 6, height = 4, dpi = 300)

message("✅ NMDS ordination by rodent species saved.")

# ------------------------------------------------------------
# Step 6: Explore Beta Diversity by Trapping Location Type
# ------------------------------------------------------------
p_nmds_location <- plot_ordination(ps_rel, nmds_ord, color = "Location_type") +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "NMDS Ordination by Trapping Location", color = "Trapping Location")

# Save NMDS plot by location
ggsave(file.path(results_dir, "Figures", "16S_Beta", "nmds_location.jpeg"),
       plot = p_nmds_location, width = 6, height = 4, dpi = 300)

message("✅ NMDS ordination by trapping location saved.")


# ************************************************************
# Script: Beta Diversity with Latitude and Longitude in 3D
# ************************************************************

message("📈 Creating 3D scatter plot for beta diversity vs latitude and longitude...")

# ------------------------------------------------------------
# Step 1: Extract NMDS Axis values
# ------------------------------------------------------------
# ************************************************************
# Script: Beta Diversity with Latitude and Longitude in 3D
# ************************************************************

message("📈 Creating 3D scatter plot for beta diversity vs latitude and longitude...")

otu_table(ps) <- otu_table(round(otu_table(ps)), taxa_are_rows = TRUE)

alpha_df <- phyloseq::estimate_richness(ps, measures = c("Observed", "Shannon", "InvSimpson")) %>%
  tibble::rownames_to_column("Sample_ID") %>%
  left_join(sample_data(ps) %>% data.frame() %>% tibble::rownames_to_column("Sample_ID"), by = "Sample_ID")

# ------------------------------------------------------------
# Step 1: Extract NMDS Axis values
# ------------------------------------------------------------
# NMDS ordination is already performed, so we extract the NMDS Axis 1 and Axis 2 values
nmds_axis1 <- nmds_ord$points[, 1]  # First NMDS axis
nmds_axis2 <- nmds_ord$points[, 2]  # Second NMDS axis

# Add NMDS axis values to the alpha_df data frame for visualization
alpha_df <- cbind(alpha_df, NMDS_axis1 = nmds_axis1, NMDS_axis2 = nmds_axis2)

# ------------------------------------------------------------
# Step 2: Create a 3D scatter plot (Latitude, Longitude, and NMDS axes)
# ------------------------------------------------------------

# Create a 3D scatter plot using plotly
p_3d <- plot_ly(alpha_df, x = ~Latitude, y = ~Longitude, z = ~NMDS_axis1,
                color = ~Morphology_species, colors = "Set3",
                type = "scatter3d", mode = "markers",
                marker = list(size = 5)) %>%
  layout(
    title = "Beta Diversity (NMDS Axis 1) vs Latitude and Longitude",
    scene = list(
      xaxis = list(title = "Latitude"),
      yaxis = list(title = "Longitude"),
      zaxis = list(title = "NMDS Axis 1 (Beta Diversity)")
    ),
    showlegend = TRUE
  )

# Display the 3D plot
p_3d

# ------------------------------------------------------------
# Step 3: Save the 3D plot as an interactive HTML file
# ------------------------------------------------------------
output_3d_file <- file.path(results_dir, "Figures", "16S_Beta", "beta_diversity_3d_latitude_longitude.html")

# Save the plot to an HTML file (interactive plot)
htmlwidgets::saveWidget(p_3d, output_3d_file)

message("✅ 3D plot saved to: ", output_3d_file)

# ------------------------------------------------------------
# Step 4: Clean up temporary objects
# ------------------------------------------------------------
rm(p_3d)

message("🧹 Cleaned up temporary objects.")

