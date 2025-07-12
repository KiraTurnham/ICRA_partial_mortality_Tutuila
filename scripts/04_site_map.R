
rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
library(purrr)
library(sf)
library(viridis)
library(cowplot)
library("svglite")

# set colors
vir_colors <- viridis(n = 4, option = "C")
print(vir_colors)

custom_colors <- vir_colors
custom_colors[4] <- "gold"  # DAA520 goldenrod 

load("data/SOUTH_COLONY_DENSITY_filtered.RData") #density in south 
load("data/ALL_COLONY_DENSITY.RData") #site data north and south unfiltered

tutuila_shape <- st_read("data/Tut_shapefiles/TUT.shp") #make sure folder has all shapefiles needed

#################################
#plot SITES across all island####
#################################
  
  #assign coordinates for PM data. Want to keep PM=0 but not NA
  ALL_site_sf <- ALL_COLONY_DENSITY %>%
  #filter(!is.na(mean_PM)) %>%
  distinct(LATITUDE, LONGITUDE, YEAR, .keep_all = TRUE) %>%  # keep all columns
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326)  # WGS84 


#plot PM
ALL_site_sf <- ALL_site_sf %>%
  mutate(shape_code = ifelse(is.na(DENSITY), "NA", as.character(YEAR)))

ALL_site_sf$YEAR <- factor(ALL_site_sf$YEAR, levels = c("2015", "2018", "2023", "2025"))

#define lat breaks from all sites_sf
lat_breaks <- ggplot_build(
  ggplot() + 
    geom_sf(data = ALL_site_sf) + 
    coord_sf()
)$layout$panel_params[[1]]$y$breaks

ggplot() +
  geom_sf(data = tutuila_shape, fill = "grey", color = "black") +
  
  # First draw all years except 2025 (background layers)
  geom_sf(data = subset(ALL_site_sf, YEAR != "2025"),
          aes(fill = YEAR, shape = YEAR),
          size = 1.5, stroke = 0.5, alpha = 0.7, color = "black") +
  
  # Then draw 2025 on top (foreground layer)
  geom_sf(data = subset(ALL_site_sf, YEAR == "2025"),
          aes(fill = YEAR, shape = YEAR),
          size = 1, stroke = 0.5, alpha = 0.7, color = "black") +
  
  scale_shape_manual(values = c("2015" = 22, "2018" = 23, "2023" = 24, "2025" = 21)) +
  scale_fill_manual(
    values = c("2015" = custom_colors[1],
               "2018" = custom_colors[2],
               "2023" = custom_colors[3],
               "2025" = custom_colors[4])
  ) +
  labs(
    x = "Longitude", y = "Latitude",
    fill = "Year",
    shape = "Year"
  ) +
  scale_y_continuous(
    breaks = lat_breaks,
    labels = ifelse(seq_along(lat_breaks) %% 2 == 1, lat_breaks, "")
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 13),# face="bold"),
    text = element_text(size = 5),
    plot.title = element_text(hjust = 0.5),
    legend.position = "top",
    legend.key.size = unit(0.4, "cm"),
    legend.box = "horizontal",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10)
  )
#ggplot2::ggsave ("plots/all_sites_by_year.png", width = 5, height = 2.5, units = 'in',  bg = "transparent")

##########################
##plot where x's are 0/NA#
##########################

# Split data##############
site_zero <- subset(ALL_site_sf, DENSITY == 0)
site_nonzero <- subset(ALL_site_sf, DENSITY != 0 | is.na(DENSITY))



ggplot() +
  geom_sf(data = tutuila_shape, fill = "grey", color = "black") +
  
  # Plot all non-2025 sites with non-zero density
  geom_sf(data = subset(site_nonzero, YEAR != "2025"),
          aes(fill = YEAR, shape = YEAR),
          size = 1.5, stroke = 0.5, alpha = 1, color = "black") +
  
  # Plot 2025 non-zero on top
  geom_sf(data = subset(site_nonzero, YEAR == "2025"),
          aes(fill = YEAR, shape = YEAR),
          size = 1, stroke = 0.5, alpha = 0.7, color = "black") +
  
  # X's for zero-density sites — drawn LAST (on top)
  geom_sf(data = site_zero,
          aes(color = YEAR),
          shape = 4, size = 1, alpha = 0.7, stroke = 0.6,
          show.legend = FALSE) +
  
  scale_shape_manual(values = c("2015" = 22, "2018" = 23, "2023" = 24, "2025" = 21)) +
  scale_fill_manual(values = c("2015" = custom_colors[1],
                               "2018" = custom_colors[2],
                               "2023" = custom_colors[3],
                               "2025" = custom_colors[4])) +
  scale_color_manual(values = c("2015" = custom_colors[1],
                                "2018" = custom_colors[2],
                                "2023" = custom_colors[3],
                                "2025" = custom_colors[4])) +
  
  labs(
    x = "Longitude", y = "Latitude",
    fill = "Year",
    shape = "Year",
    color = "Year (Density = 0)"
  ) +
  scale_y_continuous(
    breaks = lat_breaks,
    labels = ifelse(seq_along(lat_breaks) %% 2 == 1, lat_breaks, "")
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),# face="bold"),
    text = element_text(size = 5),
    plot.title = element_text(hjust = 0.5),
    legend.position = "top",
    legend.key.size = unit(0.4, "cm"),
    legend.box = "horizontal",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10)
  )



ggplot2::ggsave ("plots/all_sites_by_year_X.png", width = 5, height = 2.5, units = 'in',  bg = "transparent")
ggplot2::ggsave ("plots/all_sites_by_year_X.pdf", width = 5, height = 2.5, units = 'in',  bg = "transparent")

#################################
# plot south sites only (Fig 1a)#
#################################

#assign coordinates for density data. Want to keep PM=0 but not NA
colnames(SOUTH_COLONY_DENSITY_filtered)
sites_sf_density_south <- SOUTH_COLONY_DENSITY_filtered %>%
  distinct(LATITUDE, LONGITUDE, YEAR, .keep_all = TRUE) %>%  # keep all columns
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326)  # WGS84 

#way to plot so that 2025 data is not on top
# Make sure YEAR is a factor with the correct level order
# Set YEAR as a factor in chronological order for the legend
sites_sf_density_south$YEAR <- factor(sites_sf_density_south$YEAR, levels = c("2015", "2018", "2023", "2025"))

ggplot() +
  geom_sf(data = tutuila_shape, fill = "grey", color = "black") +
  
  # First draw all years except 2025 (background layers)
  geom_sf(data = subset(sites_sf_density_south, YEAR != "2025"),
          aes(fill = YEAR, shape = YEAR),
          size = 1.5, stroke = 0.5, alpha = 0.8, color = "black") +
  
  # Then draw 2025 on top (foreground layer)
  geom_sf(data = subset(sites_sf_density_south, YEAR == "2025"),
          aes(fill = YEAR, shape = YEAR),
          size = 2, stroke = 0.5, alpha = 0.9, color = "black") +
  
  scale_shape_manual(values = c("2015" = 22, "2018" = 23, "2023" = 24, "2025" = 21)) +
  scale_fill_manual(
    values = c("2015" = custom_colors[1],
               "2018" = custom_colors[2],
               "2023" = custom_colors[3],
               "2025" = custom_colors[4])
  ) +
  labs(
    x = "Longitude", y = "Latitude",
    fill = "Year",
    shape = "Year"
  ) +
  scale_y_continuous(
    breaks = lat_breaks,
    labels = ifelse(seq_along(lat_breaks) %% 2 == 1, lat_breaks, "")
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 14), #face="bold"),
    text = element_text(size = 6),
    plot.title = element_text(hjust = 0.5),
    legend.position = "top",
    legend.key.size = unit(0.4, "cm"),
    legend.box = "horizontal",
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 11))

ggplot2::ggsave ("paper/sites_map.pdf", width = 5, height = 3, units = 'in',  bg = "transparent")
ggplot2::ggsave ("paper/sites_map.png", width = 5, height = 3, units = 'in',  bg = "transparent")

