# ***********************************************************
# Title: Clean Rodent Field Data
# Purpose: This script cleans rodent trapping data, corrects errors,
#          standardizes date formats, and handles missing values.
# Input: "data/metadata/Rodents_catching_data.csv"
# Output: "data/processed/cleaned_rodent_data.csv"
# ***********************************************************

# Load raw rodent data
rodent_data <- read_csv(file.path(field_raw_tracking, "Rodents_catching_data.csv"))

# Fix species name inconsistencies
rodent_data$Morphology_species <- str_replace_all(rodent_data$Morphology_species,
                                                  c("Vole\\?" = "Vole", "House_mouse\\?" = "House_mouse"))

# Correct incorrect dates (if future years detected)
rodent_data$Date <- str_replace_all(rodent_data$Date, c("202[4-9]|2030" = "2023"))

# Handle missing values
rodent_data <- rodent_data %>%
  mutate(Sex = ifelse(Sex == "", "unidentified", Sex)) %>%
  mutate(Sex = gsub("63", "unidentified", Sex))


message("✅ Rodent field data cleaned successfully!")
