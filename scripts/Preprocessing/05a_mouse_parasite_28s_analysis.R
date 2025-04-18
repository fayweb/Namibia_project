# ***********************************************************
# Script: 05a_mouse_parasite_28s_analysis.R
# Purpose: Analyze 28S EMU output from Emanuel Heitlinger (rodents only)
#          and produce a parasite-focused heatmap
# ***********************************************************
#  Define file paths
counts_file <- file.path(data_dir, "processed", "28S_Sequencing_Emanuel",
                         "emu-combined-tax_id-counts.tsv")

sample_metadata_file <- file.path(data_dir, "processed",
                                  "28S_Sequencing_Emanuel", "AnimalData28S.csv")


# ***********************************************************
# Step 1: Load and filter sample metadata
# ***********************************************************
animal_28S <- read.csv(sample_metadata_file) %>%
  mutate(
    Barcode = str_trim(Barcode),
    barcode_28s = paste0("barcode", str_remove(Barcode, "^BC")),
    sample_type = str_to_lower(str_trim(sample_type)),
    project = str_to_lower(str_trim(project))
  ) %>%
  filter(
    sample_type == "rodent",
    str_starts(sample_ID, "R"),
    str_starts(project, "rodent")
  )

message("✅ Filtered to rodent samples: ", nrow(animal_28S))

# ***********************************************************
# Step 2: Reshape EMU 28S output to long format
# ***********************************************************

emu <- read.delim(counts_file)

emu_long <- emu %>%
  pivot_longer(
    cols = starts_with("barcode"),
    names_to = "barcode_28s",
    values_to = "count_28s"
  ) %>%
  filter(!is.na(count_28s) & count_28s > 0)

# ***********************************************************
# Step 3: Join with rodent metadata
# ***********************************************************

emu_annotated <- emu_long %>%
  left_join(animal_28S %>%
              dplyr::select(c("sample_ID", "sample_ID", "Date_28s_PCR",
                              "Gene", "barcode_28s")), by = "barcode_28s") %>%
  dplyr::rename(Sample_ID = sample_ID)

message("✅ Annotated EMU 28S output: ", nrow(emu_annotated), " rows")



# Join 28S annotated EMU output with rodent-level metadata
rodent_data <- rodent_data %>%
  full_join(emu_annotated, by = intersect(colnames(emu_annotated), colnames(rodent_data)))


# Check results
message("✅ Merged 28S EMU data with rodent metadata: ", nrow(rodent_data), " rows")

rm(emu_annotated, emu, emu_long, animal_28S)

