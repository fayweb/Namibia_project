# ***********************************************************
# Title: Field Nanodrop QC - Rodent Samples (March 2023)
# Purpose: Archive and assess DNA quality using Nanodrop measurements
#          collected *in the field* during rodent sampling in Namibia.
#          These were later reassessed in a laboratory in Berlin.
#
# Author: Fay Webster
# Input: Two .tsv files from Nanodrop runs (March 26–27, 2023)
# Output: Merged .csv file and exploratory quality plots
# ***********************************************************

# Create output directory for QC plots if it doesn't exist
dir.create("results/figures/archived_qc", recursive = TRUE, showWarnings = FALSE)

# ***********************************************************
# Step 1: Read and Merge Raw Nanodrop Field Files
# ***********************************************************
rodent1 <- read_tsv("data/raw/Nanodrop_measurements/Rodents_26032023.tsv")
rodent2 <- read_tsv("data/raw/Nanodrop_measurements/Rodents_27032023.tsv")

# Remove the first column (Nanodrop index column)
rodent1 <- rodent1[, -1]
rodent2 <- rodent2[, -1]

# Merge both measurement files
nanodrop_field <- bind_rows(rodent1, rodent2)

# Save combined data
write_csv(nanodrop_field, "data/processed/nanodrop_field_combined_2023.csv")
message("✅ Saved merged Nanodrop field data to: data/processed/nanodrop_field_combined_2023.csv")

# ***********************************************************
# Step 2: Clean Column Names
# ***********************************************************
# Standardize column names (replace spaces with underscores)
colnames(nanodrop_field) <- gsub(" ", "_", colnames(nanodrop_field))

# Rename specific columns for clarity
nanodrop_field <- nanodrop_field %>%
  rename(
    Quality_260_280 = `260/280`,
    Quality_260_230 = `260/230`,
    DNA_ng_ul = `Nucleic_Acid`
  )

# ***********************************************************
# Step 3: Visual QC - 260/280 Ratio vs DNA Concentration
# ***********************************************************
p1 <- ggplot(nanodrop_field, aes(x = DNA_ng_ul, y = Quality_260_280)) +
  geom_jitter() +
  labs(x = "DNA Concentration (ng/µl)", y = "260/280 Ratio",
       title = "Nanodrop: 260/280 Ratio vs DNA Concentration")
ggsave("results/figures/archived_qc/260_280_vs_concentration.png", p1, width = 8, height = 6, dpi = 300)
p1

# Filter biologically plausible 260/280 values
p2 <- nanodrop_field %>%
  filter(Quality_260_280 < 2.5, Quality_260_280 > 1.5) %>%
  ggplot(aes(x = DNA_ng_ul, y = Quality_260_280)) +
  geom_jitter() +
  labs(x = "DNA Concentration (ng/µl)", y = "260/280 Ratio",
       title = "Filtered: 260/280 Ratio (1.5–2.5)")
ggsave("results/figures/archived_qc/260_280_filtered.png", p2, width = 8, height = 6, dpi = 300)
p2

# ***********************************************************
# Step 4: Visual QC - 260/230 Ratios
# ***********************************************************
p3 <- ggplot(nanodrop_field, aes(x = DNA_ng_ul, y = Quality_260_230)) +
  geom_jitter() +
  labs(x = "DNA Concentration (ng/µl)", y = "260/230 Ratio",
       title = "Nanodrop: 260/230 Ratio vs DNA Concentration")
ggsave("results/figures/archived_qc/260_230_vs_concentration.png", p3, width = 8, height = 6, dpi = 300)
p3

# Filter out extreme values (< 60)
p4 <- nanodrop_field %>%
  filter(Quality_260_230 < 60) %>%
  ggplot(aes(x = DNA_ng_ul, y = Quality_260_230)) +
  geom_jitter() +
  labs(x = "DNA Concentration (ng/µl)", y = "260/230 Ratio",
       title = "Filtered: 260/230 Ratio (< 60)")
ggsave("results/figures/archived_qc/260_230_filtered.png", p4, width = 8, height = 6, dpi = 300)
p4

# Focus on optimal purity range (1.8–2.5)
p5 <- nanodrop_field %>%
  filter(Quality_260_230 < 2.5, Quality_260_230 > 1.8) %>%
  ggplot(aes(x = DNA_ng_ul, y = Quality_260_230)) +
  geom_jitter() +
  labs(x = "DNA Concentration (ng/µl)", y = "260/230 Ratio",
       title = "Filtered: 260/230 Ratio (1.8–2.5)")
ggsave("results/figures/archived_qc/260_230_ideal_range.png", p5, width = 8, height = 6, dpi = 300)
p5

# ***********************************************************
# Step 5: Filter “Golden Ratio” DNA Samples
# ***********************************************************
# These samples had clean 260/280 ratios (indicative of good DNA)
golden <- nanodrop_field %>%
  filter(Quality_260_280 < 2.5, Quality_260_280 > 1.5)

write_csv(golden, "data/processed/nanodrop_field_golden_2023.csv")
message("✅ Filtered golden ratio samples saved to: data/processed/nanodrop_field_golden_2023.csv")

# ***********************************************************
# Notes:
# - These measurements were taken directly in the field during March 2023.
# - DNA quality was later reassessed in the lab in Berlin under controlled conditions.
# - These field values are preserved here for provenance and exploratory comparison only.
# ***********************************************************

