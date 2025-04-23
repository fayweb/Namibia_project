# ************************************************************
# Script: 6_heatmap.R
# Purpose: Create heatmaps for taxa prevalence and abundance
# Input:   processed/phyloseq_16s_decontaminated.rds
# Output:  Results/Figures/16S_QC/heatmap_plots.jpeg
# ************************************************************

message("🔬 Starting heatmap creation...")

# ------------------------------------------------------------
# Step 1: Load the phyloseq object and apply filters
# ------------------------------------------------------------
ps <- readRDS(file.path(processed_data, "phyloseq_16s_cleaned_genera.rds"))

# ------------------------------------------------------------
# Step 2: Apply tax_fix() to standardize the taxonomy
# ------------------------------------------------------------
ps <- ps %>%
  tax_fix(min_length = 4, unknowns = c(""), sep = " ", anon_unique = TRUE, suffix_rank = "classified") %>%
  phyloseq_validate()

# ------------------------------------------------------------
# Step 3: Taxonomy Cleaning and Preprocessing for Heatmap
# ------------------------------------------------------------
# Copy of the phyloseq object for manipulation
psq <- ps

# Filter taxa by minimum prevalence
psq <- tax_filter(psq, min_prevalence = 5)

# Aggregate taxa by genus
psq <- tax_agg(psq, rank = "genus")

# Get top 20 genera based on prevalence
taxa <- tax_top(psq, n = 20, rank = "genus")

# ------------------------------------------------------------
# Step 4: Prevalence Barplot Annotation
# ------------------------------------------------------------
heatmapAnnoFunction <- anno_tax_prev()(data = psq, which = "row", taxa = taxa)

# Draw the prevalence barplot
vp <- viewport(width = 0.75, height = 0.75)
grid::grid.newpage()
pushViewport(vp)
draw(heatmapAnnoFunction)

# ------------------------------------------------------------
# Step 5: Abundance Boxplot Annotation
# ------------------------------------------------------------
#psq <- tax_fix(psq)  # Ensure that taxonomic issues are fixed
heatmapAnnoFunction <- anno_tax_box()(data = psq, which = "column", taxa = taxa)

# Draw the abundance boxplot
vp <- viewport(width = 0.75, height = 0.75)
grid.newpage()
pushViewport(vp)
draw(heatmapAnnoFunction)

# ------------------------------------------------------------
# Step 6: Generate the Final Heatmap with Top Taxa
# ------------------------------------------------------------
# Transform the data into compositional format
psq2 <- tax_transform(psq, "compositional", rank = "genus")

# Create the heatmap for the top 20 taxa
best_heatmap <- comp_heatmap(
  data = psq2,
  taxa = taxa[1:20],
  tax_anno = taxAnnotation(Prev = anno_tax_prev(undetected = 50))
)

# Draw the heatmap
best_heatmap_2 <- ComplexHeatmap::draw(object = best_heatmap)

# Draw the heatmap and save it to a file
png(file.path(results_dir, "Figures", "16S_QC", "heatmap_plots.png"), width = 8, height = 5, units = "in", res = 300)
ComplexHeatmap::draw(best_heatmap)
dev.off()  # Close the graphics device


message("✅ Heatmap plot saved.")


# ------------------------------------------------------------
# Final Step: Clean up temporary objects
# ------------------------------------------------------------
rm(list = c("ps", "psq", "taxa", "best_heatmap", "best_heatmap_2"))

message("🧹 Cleaned up temporary objects.")

