library(dplyr)
library(tidyr)
library(ggplot2)
library(tmap)
library(tmaptools)
library(sf)
library(rnaturalearth)
library(sf)
library(cowplot)
library(magick)


df <- read.csv("Data/Field_tracking_data/Rodents_catching_data.csv")

glimpse(df)

# correct mistakes
df$Morphology_species <- gsub("Vole\\?", "Vole", df$Morphology_species)

# convert the df into a spatial df
df_sf <- st_as_sf(df, coords = c("Longitude", "Latitude"), crs = 4326)

#tmap_mode("view")  # Use this for interactive maps

# map on species and geolocations
tm_basemap("OpenStreetMap") + 
    tm_shape(df_sf) + 
    tm_dots(col = "Morphology_species", 
            title = "Species", 
            palette = "Set1", 
            size = 0.5, 
            border.col = "white", 
            id = "Sample_ID")

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
    rel_widths = c(1, 0.1, 1) # Relative widths of the columns. Adjust the 0.1 to increase/decrease spacing.
)

print(combined_plot)


ggsave("Lodge_location.png", combined_plot, width = 10, height = 7, dpi = 300)


# Location: Bildah
img1 <- ggdraw() + draw_image("Figures/Maps/3.Bildah_species.png")
img2 <- ggdraw() + draw_image("Figures/Maps/4.Bildah_type.png")

combined_plot <- plot_grid(
    img1, NULL, img2,         # The NULL creates an empty column for spacing
    labels = c("C. Rodent species", "", "D. Location of rodents"), # Labels for each panel. The empty label ("") corresponds to the NULL/spacing column.
    label_size = 14,          # Adjust size as needed
    ncol = 3,                 # Increase to 3 because of the added NULL/spacing column
    rel_widths = c(1, 0.1, 1) # Relative widths of the columns. Adjust the 0.1 to increase/decrease spacing.
)

print(combined_plot)



ggsave("Bildah_locations.png", combined_plot, width = 10, height = 7, dpi = 300)

# Location: Cheetah research station
img1 <- ggdraw() + draw_image("Figures/Maps/3.Cheetah_research_station.png")
img2 <- ggdraw() + draw_image("Figures/Maps/6.Cheetah_research_station_type.png")

combined_plot <- plot_grid(
    img1, NULL, img2,         # The NULL creates an empty column for spacing
    labels = c("E. Rodent species", "", "F. Location of rodents"), # Labels for each panel. The empty label ("") corresponds to the NULL/spacing column.
    label_size = 14,          # Adjust size as needed
    ncol = 3,                 # Increase to 3 because of the added NULL/spacing column
    rel_widths = c(1, 0.1, 1) # Relative widths of the columns. Adjust the 0.1 to increase/decrease spacing.
)

print(combined_plot)


ggsave("Cheetah_location.png", combined_plot, width = 10, height = 7, dpi = 300)

