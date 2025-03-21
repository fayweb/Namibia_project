# ***********************************************************
# Title: Composite Site Panel Maps
# Purpose: Combine map figures into labeled panels
# Output: Saves final PNG panels in results/figures/maps
# ***********************************************************

# Okambara
img1 <- ggdraw() + draw_image("results/figures/maps/2.Elephant_lodge_location_species.png")
img2 <- ggdraw() + draw_image("results/figures/maps/5.Elephant_lodge_location_type.png")

combined_plot <- plot_grid(
  img1, NULL, img2,
  labels = c("A. Rodent species", "", "B. Location of rodents"),
  label_size = 14, ncol = 3, rel_widths = c(1, 0.1, 1), label_y = 0.8
)
ggsave("results/figures/maps/Lodge_location.png", combined_plot, width = 10, height = 7, dpi = 300)

# Bildah
img1 <- ggdraw() + draw_image("results/figures/maps/3.Bildah_species.png")
img2 <- ggdraw() + draw_image("results/figures/maps/4.Bildah_type.png")
combined_plot <- plot_grid(img1, NULL, img2,
                           labels = c("C. Rodent species", "", "D. Location of rodents"),
                           label_size = 14, ncol = 3, rel_widths = c(1, 0.1, 1), label_y = 0.8
)
ggsave("results/figures/maps/Bildah_locations.png", combined_plot, width = 10, height = 7, dpi = 300)

# Cheetah Station
img1 <- ggdraw() + draw_image("results/figures/maps/3.Cheetah_research_station.png")
img2 <- ggdraw() + draw_image("results/figures/maps/6.Cheetah_research_station_type.png")
combined_plot <- plot_grid(img1, NULL, img2,
                           labels = c("E. Rodent species", "", "F. Location of rodents"),
                           label_size = 14, ncol = 3, rel_widths = c(1, 0.1, 1), label_y = 0.8
)
ggsave("results/figures/maps/Cheetah_location.png", combined_plot, width = 10, height = 7, dpi = 300)

message("✅ Site panel figures saved!")

