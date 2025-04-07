# Load required libraries
library(tidyverse)
library(ggplot2)
library(janitor)

# Read merged OTU + taxonomy + metadata files
otu_marly <- read.csv("Data/processed/EMU_output/marly_standard_filtering/otu_taxonomy_metadata_marly.csv")
otu_melanie <- read.csv("Data/processed/EMU_output/melanie_lenient_filtering/otu_taxonomy_metadata_melanie.csv")

# Quick check of structure
glimpse(otu_marly)
glimpse(otu_melanie)


# Add a filtering label to each dataset
otu_marly$filtering <- "Marly_standard"
otu_melanie$filtering <- "Melanie_lenient"

# Bind both datasets
otu_combined <- dplyr::bind_rows(otu_marly, otu_melanie)

# Summarize total reads per sample
read_depth <- otu_combined %>%
  dplyr::group_by(barcode, filtering) %>%
  dplyr::summarise(total_reads = sum(count), .groups = "drop")

# Plot read depth
ggplot(read_depth, aes(x = reorder(barcode, -total_reads), y = total_reads, fill = filtering)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Read Depth per Sample by Filtering Strategy",
       x = "Sample Barcode",
       y = "Total Reads") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
