library(dplyr)
library(tidyr)
library(ggplot2)
library(mapview)
library(sf)
library(rnaturalearth)
library(sf)
library(cowplot)
library(magick)
library(RColorBrewer)
library(stringr)
library(tmap)

# set the current repository to project directory before
#setwd("~/GitHub/Namibia_project/")

df <- read.csv("Data/Field_tracking_data/Rodents_catching_data.csv")

glimpse(df)

# correct mistakes
df$Morphology_species <- 
    gsub("Vole\\?", "Vole", df$Morphology_species)

df$Morphology_species <- 
    gsub("House_mouse\\?", "House_mouse", df$Morphology_species)

df$Date <- str_replace(df$Date, "2024", "2023")
df$Date <- str_replace(df$Date, "2025", "2023")
df$Date <- str_replace(df$Date, "2026", "2023")
df$Date <- str_replace(df$Date, "2027", "2023")
df$Date <- str_replace(df$Date, "2028", "2023")
df$Date <- str_replace(df$Date, "2029", "2023")
df$Date <- str_replace(df$Date, "2030", "2023")


df <- df %>%
    mutate(Sex = ifelse(Sex == "", "unidentified", Sex))

df$Sex <- 
    gsub("63", "unidentified", df$Sex)

# convert the df into a spatial df
df_sf <- st_as_sf(df, coords = c("Longitude", "Latitude"), crs = 4326)


unique_species <- unique(df_sf$Morphology_species)
color_scheme <- brewer.pal(length(unique_species), "Set1")
names(color_scheme) <- unique_species


# map on species and geolocations
mapview(df_sf, 
        zcol = "Morphology_species", 
        col.regions = color_scheme, 
        map.types = "OpenStreetMap",
        popup = "Sample_ID",
        legend = TRUE)


# map on location and structure where the rodent was caught
tm_basemap("OpenStreetMap") + 
    tm_shape(df_sf) + 
    tm_dots(col = "Location_type", 
            title = "Species", 
            palette = "Set1", 
            size = 0.5, 
            border.col = "white", 
            id = "Sample_ID")


# location: Okambara Elephan Lodge
img1 <- ggdraw() + draw_image("Figures/Maps/2.Elephant_lodge_location_species.png")
img2 <- ggdraw() + draw_image("Figures/Maps/5.Elephant_lodge_location_type.png")

combined_plot <- plot_grid(
    img1, NULL, img2,         # The NULL creates an empty column for spacing
    labels = c("A. Rodent species", "", "B. Location of rodents"), # Labels for each panel. The empty label ("") corresponds to the NULL/spacing column.
    label_size = 14,          # Adjust size as needed
    ncol = 3,                 # Increase to 3 because of the added NULL/spacing column
    rel_widths = c(1, 0.1, 1), # Relative widths of the columns. Adjust the 0.1 to increase/decrease spacing.
    label_y = 0.8       # Adjust this value to move the label closer to the plot
)

print(combined_plot)


ggsave("Figures/Maps/Lodge_location.png", combined_plot, width = 10, height = 7, dpi = 300)


# Location: Bildah
img1 <- ggdraw() + draw_image("Figures/Maps/3.Bildah_species.png")
img2 <- ggdraw() + draw_image("Figures/Maps/4.Bildah_type.png")

combined_plot <- plot_grid(
    img1, NULL, img2,         # The NULL creates an empty column for spacing
    labels = c("C. Rodent species", "", "D. Location of rodents"), # Labels for each panel. The empty label ("") corresponds to the NULL/spacing column.
    label_size = 14,          # Adjust size as needed
    ncol = 3,                 # Increase to 3 because of the added NULL/spacing column
    rel_widths = c(1, 0.1, 1), # Relative widths of the columns. Adjust the 0.1 to increase/decrease spacing.
    label_y = 0.8)

print(combined_plot)



ggsave("Figures/Maps/Bildah_locations.png", combined_plot, width = 10, height = 7, dpi = 300)

# Location: Cheetah research station
img1 <- ggdraw() + draw_image("Figures/Maps/3.Cheetah_research_station.png")
img2 <- ggdraw() + draw_image("Figures/Maps/6.Cheetah_research_station_type.png")

combined_plot <- plot_grid(
    img1, NULL, img2,         # The NULL creates an empty column for spacing
    labels = c("E. Rodent species", "", "F. Location of rodents"), # Labels for each panel. The empty label ("") corresponds to the NULL/spacing column.
    label_size = 14,          # Adjust size as needed
    ncol = 3,                 # Increase to 3 because of the added NULL/spacing column
    rel_widths = c(1, 0.1, 1), # Relative widths of the columns. Adjust the 0.1 to increase/decrease spacing.
    label_y = 0.8)

print(combined_plot)


ggsave("Figures/Maps/Cheetah_location.png", combined_plot, width = 10, height = 7, dpi = 300)


# How many males?
df %>%
    group_by(Sex) %>% 
    summarise(n())

# Average weight per species
figure_df <- df %>%
    group_by(Morphology_species) %>%
    summarise(mean_weight = mean(Weight_g, na.rm = TRUE)) %>%
    arrange(-mean_weight)  # Arrange by mean_weight in descending order


# To view the new data frame
View(figure_df)

figure_df$Morphology_species <- factor(figure_df$Morphology_species, levels = figure_df$Morphology_species)

# Plot
p <- ggplot(figure_df, aes(x = Morphology_species, y = mean_weight)) +
    geom_bar(stat="identity", fill="orchid1") +
    labs(title = "Average Weight per Species", x = "Species", y = "Average Weight (g)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot as an image
ggsave("Figures/Plots/average_weight_per_species.png", p, width = 10, height = 7, dpi = 300)


# Sex counts
sex_counts <- df %>%
    filter(!is.na(Sex)) %>% # Filter out any rows with missing Sex data
    group_by(Sex) %>%
    summarise(count = n())
# Plot
p <- ggplot(sex_counts, aes(x = "", y = count, fill = Sex)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0) + # Transform the bar chart into a pie chart
    labs(title = "Distribution of Mice by Sex", fill = "Sex") +
    theme_minimal() +
    theme(axis.text.x = element_blank(), # Hide x-axis labels
          axis.title.x = element_blank(), # Hide x-axis title
          axis.title.y = element_blank()) # Hide y-axis title

# Save the plot as an image
ggsave("Figures/Plots/sex_distribution_pie_chart.png", p, width = 7, height = 7, dpi = 300)


