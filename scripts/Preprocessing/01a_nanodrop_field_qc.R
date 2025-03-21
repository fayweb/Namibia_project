# ***********************************************************
# Title: Field Nanodrop QC - Rodent Samples (March 2023)
# Purpose: Read, clean, and assess Nanodrop measurements
# Authors: Fay Webster
# Input: Two .tsv files from Nanodrop runs
# Output: Merged .csv and QC plots
# ***********************************************************

# Create output directory for QC plots if it doesn't exist
dir.create("results/figures/archived_qc", recursive = TRUE, showWarnings = FALSE)

# ***********************************************************
# Step 1: Read and Merge Nanodrop Files
# ***********************************************************
rodent1 <- read_tsv("data/raw/Nanodrop_measurements/Rodents_26032023.tsv")
rodent2 <- read_tsv("data/raw/Nanodrop_measurements/Rodents_27032023.tsv")

# Remove first column (ID column from Nanodrop export)
rodent1 <- rodent1[, -1]
rodent2 <- rodent2[, -1]

# Combine
nanodrop_field <- bind_rows(rodent1, rodent2)

# Save merged file
write_csv(nanodrop_field, "data/processed/Nanodrop_rodents_202303_combined.csv")
message("✅ Saved merged Nanodrop field data to: data/processed/Nanodrop_rodents_202303_combined.csv")

# ***********************************************************
# Step 2: Tidy Column Names (Optional but recommended)
# ***********************************************************
colnames(nanodrop_field) <- gsub(" ", "_", colnames(nanodrop_field))

# Rename columns for clarity if needed
nanodrop_field <- nanodrop_field %>%
  rename(
    Quality_260_280 = `260/280`,
    Quality_260_230 = `260/230`,
    DNA_ng_ul = `Nucleic_Acid`
  )

# ***********************************************************
# Step 3: Plot 260/280 Ratio vs DNA Concentration
# ***********************************************************
p1 <- ggplot(nanodrop_field, aes(x = DNA_ng_ul, y = Quality_260_280)) +
  geom_jitter() +
  labs(x = "DNA Concentration (ng/µl)", y = "260/280 Ratio",
       title = "Nanodrop: 260/280 Ratio vs DNA Concentration")
ggsave("results/figures/archived_qc/260_280_vs_concentration.png", p1, width = 8, height = 6, dpi = 300)

p1
# Filter outliers (extreme 260/280)
p2 <- nanodrop_field %>%
  filter(Quality_260_280 < 2.5, Quality_260_280 > 1.5) %>%
  ggplot(aes(x = DNA_ng_ul, y = Quality_260_280)) +
  geom_jitter() +
  labs(x = "DNA Concentration (ng/µl)", y = "260/280 Ratio",
       title = "Filtered: 260/280 Ratio (1.5 - 2.5)")
ggsave("results/figures/archived_qc/260_280_filtered.png", p2, width = 8, height = 6, dpi = 300)


p2
# ***********************************************************
# Step 4: Plot 260/230 Ratios
# ***********************************************************
p3 <- ggplot(nanodrop_field, aes(x = DNA_ng_ul, y = Quality_260_230)) +
  geom_jitter() +
  labs(x = "DNA Concentration (ng/µl)", y = "260/230 Ratio",
       title = "Nanodrop: 260/230 Ratio vs DNA Concentration")
ggsave("results/figures/archived_qc/260_230_vs_concentration.png", p3, width = 8, height = 6, dpi = 300)

p3

# Filter high 260/230 values
p4 <- nanodrop_field %>%
  filter(Quality_260_230 < 60) %>%
  ggplot(aes(x = DNA_ng_ul, y = Quality_260_230)) +
  geom_jitter() +
  labs(x = "DNA Concentration (ng/µl)", y = "260/230 Ratio",
       title = "Filtered: 260/230 Ratio (<60)")
ggsave("results/figures/archived_qc/260_230_filtered.png", p4, width = 8, height = 6, dpi = 300)

p4
# Focus on ideal 260/230 range
p5 <- nanodrop_field %>%
  filter(Quality_260_230 < 2.5, Quality_260_230 > 1.8) %>%
  ggplot(aes(x = DNA_ng_ul, y = Quality_260_230)) +
  geom_jitter() +
  labs(x = "DNA Concentration (ng/µl)", y = "260/230 Ratio",
       title = "Filtered: 260/230 Ratio (1.8 - 2.5)")
ggsave("results/figures/archived_qc/260_230_ideal_range.png", p5, width = 8, height = 6, dpi = 300)

p5
# ***********************************************************
# Step 5: Filter Golden Ratio Samples
# ***********************************************************
golden <- nanodrop_field %>%
  filter(Quality_260_280 < 2.5, Quality_260_280 > 1.5)

write_csv(golden, "data/processed/Nanodrop_rodents_202303_golden_ratio.csv")
message("✅ Filtered golden ratio samples saved to: data/processed/Nanodrop_rodents_202303_golden_ratio.csv")

