# ***********************************************************
# Title: Temporal Patterns Analysis
# Purpose: Analyzes species occurrences over time
# Input: "data/processed/cleaned_rodent_data.csv"
# Output: "results/figures/temporal_patterns.png"
# ***********************************************************

# Load cleaned dataset
rodent_data <- read_csv("data/processed/cleaned_rodent_data.csv")

# Convert Date to date format
rodent_data$Date <- dmy(rodent_data$Date)

# Generate histogram of species occurrences over time
p <- ggplot(rodent_data, aes(x = Date, fill = Morphology_species)) +
  geom_histogram(binwidth = 1) +
  labs(title = "Sampling Counts by Date", x = "Date", y = "Count") +
  theme_minimal()

p

# Save plot
ggsave("results/figures/temporal_patterns.png", p, width = 10, height = 7, dpi = 300)
message("✅ Temporal patterns analysis saved!")
