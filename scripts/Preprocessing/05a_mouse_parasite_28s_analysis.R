# ***********************************************************
# Script: 05a_mouse_parasite_28s_analysis.R
# Purpose: Analyze 28S EMU output from Emanuel Heitlinger (rodents only)
#          and produce a parasite-focused heatmap
# ***********************************************************

# 📂 Define file paths
counts_file <- file.path(data_dir, "emu-combined-tax_id-counts.tsv")
sample_metadata_file <- file.path(data_dir, "AnimalData28S.csv")
output_fig <- file.path("Results", "Figures", "rodents_parasite_heatmap.png")

write.csv(x = emu, file = "Data/processed/28S_Sequencing_Emanuel/emu_combined_tx_id_counts.csv", row.names = FALSE)
#  Read EMU combined abundance counts
emu <- read.delim("Data/processed/28S_Sequencing_Emanuel/emu-combined-tax_id-counts.tsv")

# Separate taxonomic data and count matrix
emuCounts <- emu %>%
  select(starts_with("barcode")) %>%
  round() %>%
  replace(is.na(.), 0)

emuTax <- emu %>%
  select(-starts_with("barcode")) %>%
  mutate_all(as.character) %>%
  mutate_all(~gsub("^$", NA, .)) %>%
  dplyr::select(tax_id, superkingdom, phylum, class, order, family, genus, species)

# 📋 Load sample metadata
samples <- read.csv("Data/processed/28S_Sequencing_Emanuel/AnimalData28S.csv")
samples$Barcode <- gsub("BC(\\d\\d) *", "barcode\\1", samples$Barcode)
rownames(samples) <- samples$Barcode

# 📌 Subset to rodent samples
samples <- samples[samples$project %in% c("rodents", "rodeants"), ]
emuCounts <- emuCounts[, colnames(emuCounts) %in% rownames(samples)]
samples <- samples[colnames(emuCounts), ]  # match order

# 🧪 Build phyloseq object
ps <- phyloseq(
  otu_table(emuCounts, taxa_are_rows = TRUE),
  tax_table(as.matrix(emuTax)),
  sample_data(samples)
)

# 🎯 Subset to relevant parasite taxa
ps_parasites <- subset_taxa(ps, phylum %in% c("Nematoda", "Apicomplexa", "Platyhelminthes"))
ps_parasites <- subset_taxa(ps_parasites, !order %in% "Eugregarinorida")  # remove Gregarines
ps_parasites <- subset_taxa(ps_parasites, taxa_sums(ps_parasites) > 1)

# 📊 Generate heatmap
message("🔹 Plotting rodent parasite heatmap...")
png(output_fig, width = 1200, height = 600)
pheatmap(
  log10(otu_table(ps_parasites) + 1),
  labels_row = tax_table(ps_parasites)[, "species"],
  labels_col = sample_data(ps_parasites)$sample_ID,
  display_numbers = FALSE,
  fontsize_row = 6,
  fontsize_col = 8
)
dev.off()

message("✅ Mouse parasite 28S heatmap saved to:", output_fig)
