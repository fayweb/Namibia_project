# ***********************************************************
# Title: Species Distribution Analysis
# Purpose: This script counts and visualizes species distribution
# Input: "data/processed/cleaned_rodent_data.csv"
# Output: "results/figures/species_distribution.png"
# ***********************************************************

# Count species occurrences
species_count <- rodent_data %>%
  group_by(Morphology_species) %>%
  summarize(count = n())

# Generate bar plot
p <- ggplot(species_count, aes(x = Morphology_species, y = count, fill = Morphology_species)) +
  geom_bar(stat = "identity") +
  labs(title = "Distribution of Rodent Species", x = "Species", y = "Count") +
  theme_minimal()

#p

# Save plot
ggsave("results/figures/species_distribution.png", p, width = 10, height = 7, dpi = 300)
message("✅ Species distribution plot saved!")
