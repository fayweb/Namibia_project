# ***********************************************************
# Title: Weight Distribution Analysis
# Purpose: Analyzes weight distribution across rodent species
# Input: "data/processed/cleaned_rodent_data.csv"
# Output: "results/figures/weight_distribution.png"
# ***********************************************************

# Load cleaned dataset
rodent_data <- read_csv("data/processed/cleaned_rodent_data.csv")

# Remove NA weights for visualization
rodent_data_clean <- rodent_data %>% drop_na(Weight_g)

# Generate violin plot for weight distribution
p <- ggplot(rodent_data_clean, aes(x = Morphology_species, y = Weight_g, fill = Morphology_species)) +
  geom_violin() +
  labs(title = "Weight Distribution by Species", x = "Species", y = "Weight (g)") +
  theme_minimal()

p



# Save plot
ggsave("results/figures/weight_distribution.png", p, width = 10, height = 7,
       dpi = 300)

# Average weight bar plot
figure_df <- rodent_data %>%
  group_by(Morphology_species) %>%
  summarise(mean_weight = mean(Weight_g, na.rm = TRUE)) %>%
  arrange(-mean_weight)

figure_df$Morphology_species <- factor(figure_df$Morphology_species, levels = figure_df$Morphology_species)

p <- ggplot(figure_df, aes(x = Morphology_species, y = mean_weight)) +
  geom_bar(stat = "identity", fill = "orchid1") +
  labs(title = "Average Weight per Species", x = "Species", y = "Average Weight (g)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p

ggsave("results/figures/average_weight_per_species.png", p, width = 10, height = 7, dpi = 300)
message("✅ Average weight plot saved!")

message("✅ Weight distribution analysis saved!")

