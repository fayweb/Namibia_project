# ***********************************************************
# Title: Habitat Distribution Analysis
# Purpose: Analyzes habitat preference of rodent species
# Input: "data/processed/cleaned_rodent_data.csv"
# Output: "results/figures/habitat_distribution.png"
# ***********************************************************

rodent_data <- read_csv("data/processed/cleaned_rodent_data.csv")

habitat_distribution <- rodent_data %>%
  group_by(Morphology_species, Location_type) %>%
  summarize(count = n(), .groups = "drop")

p <- ggplot(habitat_distribution, aes(x = Location_type, y = count, fill = Morphology_species)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Species Distribution by Habitat Type", x = "Habitat Type", y = "Count") +
  theme_minimal()

ggsave("results/figures/habitat_distribution.png", p, width = 10, height = 7, dpi = 300)
message("✅ Habitat distribution analysis saved!")
