# ***********************************************************
# Title: Sex & Age Distribution
# Purpose: Analyzes and visualizes sex & age distribution
# Input: "data/processed/cleaned_rodent_data.csv"
# Output: "results/figures/sex_age_distribution.png"
# ***********************************************************
# Load cleaned dataset
rodent_data <- read_csv("data/processed/cleaned_rodent_data.csv")

# Group by species, sex, and age
sex_age_distribution <- rodent_data %>%
  group_by(Morphology_species, Sex, Age) %>%
  summarize(count = n()) %>%
  ungroup()

# Generate bar plot
p <- ggplot(sex_age_distribution, aes(x = Morphology_species, y = count, fill = Sex)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Sex Distribution by Species", x = "Species", y = "Count") +
  theme_minimal()

p

ggsave("results/figures/sex_age_distribution.png", p, width = 10, height = 7,
       dpi = 300)

# Sex pie chart
sex_counts <- rodent_data %>%
  filter(!is.na(Sex)) %>%
  group_by(Sex) %>%
  summarise(count = n())

p_pie <- ggplot(sex_counts, aes(x = "", y = count, fill = Sex)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  labs(title = "Distribution of Mice by Sex", fill = "Sex") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


p_pie

ggsave("results/figures/sex_distribution_pie_chart.png", p_pie, width = 7, height = 7, dpi = 300)
message("✅ Sex pie chart saved!")

