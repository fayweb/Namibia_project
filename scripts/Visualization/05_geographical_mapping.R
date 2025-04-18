# ***********************************************************
# Title: Species Geographical Distribution
# Purpose: Maps species locations using Leaflet
# Input: "data/processed/cleaned_rodent_data.csv"
# Output: Interactive Map (Leaflet)
# ***********************************************************

# ***********************************************************
# Step 1: Load & Prepare Data
# ***********************************************************


# Ensure Longitude and Latitude are numeric
rodent_data <- rodent_data %>%
  mutate(Longitude = as.numeric(Longitude),
         Latitude = as.numeric(Latitude)) %>%
  filter(!is.na(Longitude) & !is.na(Latitude))  # Remove missing coordinates

# ***********************************************************
# Step 2: Convert Data to Spatial Format (sf)
# ***********************************************************

df_sf <- st_as_sf(rodent_data, coords = c("Longitude", "Latitude"), crs = 4326)

# Extract Longitude & Latitude back into separate columns for Leaflet compatibility
df_sf <- df_sf %>%
  mutate(Longitude = st_coordinates(.)[, 1],
         Latitude = st_coordinates(.)[, 2])

# ***********************************************************
# Step 3: Generate Interactive Map with Leaflet
# ***********************************************************

m <- leaflet(df_sf) %>%
  addTiles() %>%
  addCircleMarkers(lng = ~Longitude, lat = ~Latitude, popup = ~Morphology_species,
                   color = "blue", radius = 3, opacity = 0.8) %>%
  addLegend("bottomright", pal = colorFactor(topo.colors(10), df_sf$Morphology_species),
            values = df_sf$Morphology_species, title = "Species")

message("✅ Geographical mapping completed!")
m  # Display the map
