# ************************************************************
# Script: 02_decontam_16s.R
# Purpose: Identify and remove contaminant taxa using `decontam`
# Input:   processed/phyloseq_16s_final.rds
# Output:  processed/phyloseq_16s_decontaminated.rds
# ************************************************************

# ------------------------------------------------------------
# Step 1: Load phyloseq object
# ------------------------------------------------------------
ps <- readRDS(file.path(processed_data, "phyloseq_16s_final.rds"))
message("✅ Loaded phyloseq object: ", nsamples(ps), " samples")

# ------------------------------------------------------------
# Step 2: Visualize library size per sample
# ------------------------------------------------------------
sample_data(ps)$LibrarySize <- sample_sums(ps)

ggplot(as.data.frame(sample_data(ps)),
       aes(x = reorder(rownames(sample_data(ps)), LibrarySize),
           y = LibrarySize,
           color = sample_data(ps)$sample_or_control)) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  labs(x = "Sample Index", y = "Library Size (reads)", title = "Read Depth by Sample Type")

# ------------------------------------------------------------
# Step 3: Identify contaminants using DNA concentration
# ------------------------------------------------------------
sample_data(ps)$is.neg <- sample_data(ps)$sample_or_control == "control"

contam_df <- isContaminant(ps,
                           method = "combined",
                           neg = "is.neg",
                           conc = "conc_16s__PCR",
                           threshold = 0.2)

message("🧪 Identified contaminants: ", sum(contam_df$contaminant), "/", nrow(contam_df))

# ------------------------------------------------------------
# Step 4: Visualize prevalence in controls vs. samples
# ------------------------------------------------------------
ps.pa <- transform_sample_counts(ps, function(x) 1*(x > 0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$sample_or_control == "control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$sample_or_control == "sample", ps.pa)

df.pa <- data.frame(
  pa.neg = taxa_sums(ps.pa.neg),
  pa.pos = taxa_sums(ps.pa.pos),
  contaminant = contam_df$contaminant
)

ggplot(df.pa, aes(x = pa.neg, y = pa.pos, color = contaminant)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  labs(title = "Prevalence: Contaminants vs. Non-Contaminants",
       x = "Prevalence in Neg. Controls",
       y = "Prevalence in True Samples")

# ------------------------------------------------------------
# Step 5: Prune contaminants from phyloseq object
# ------------------------------------------------------------
contaminant_taxa <- rownames(contam_df)[contam_df$contaminant]
ps_decont <- prune_taxa(!taxa_names(ps) %in% contaminant_taxa, ps)

message("🧼 Removed contaminants — new taxa count: ", ntaxa(ps_decont))

# ------------------------------------------------------------
# Step 6: Save cleaned phyloseq object
# ------------------------------------------------------------
decont_output_path <- file.path(processed_data, "phyloseq_16s_decontaminated.rds")
saveRDS(ps_decont, file = decont_output_path)
message("💾 Saved decontaminated phyloseq object to: ", decont_output_path)

# ------------------------------------------------------------
# Step 7: Clean up
# ------------------------------------------------------------
rm(contam_df, df.pa, ps.pa, ps.pa.neg, ps.pa.pos, contaminant_taxa)
message("🧹 Cleaned up temporary objects.")
