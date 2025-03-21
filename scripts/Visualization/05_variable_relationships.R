# ***********************************************************
# Title: Variable Relationships (Geographical Scatter Plot)
# Purpose: Analyzes relationships between species and geographical location
# Input: "data/processed/cleaned_rodent_data.csv"
# Output: "results/figures/species_scatter_plot.png"
# ***********************************************************
# Load cleaned dataset
rodent_data <- read_csv("data/processed/cleaned_rodent_data.csv")

# Generate scatter plot of species distribution across longitude & latitude
p <- ggplot(rodent_data, aes(x = Longitude, y = Latitude, color = Morphology_species)) +
  geom_point() +
  labs(title = "Species Distribution by Geographic Location", x = "Longitude", y = "Latitude") +
  theme_minimal()
p

# Save plot
ggsave("results/figures/species_scatter_plot.png", p, width = 10, height = 7, dpi = 300)
message("✅ Variable relationships (scatter plot) saved!")
