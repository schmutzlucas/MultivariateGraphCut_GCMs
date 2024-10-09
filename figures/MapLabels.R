# Generate the polychrome color palette with 26 colors
color_palette <- pals::glasbey(length(model_names))

# Correct latitudes: adjust for the fact that your latitude values are from -70 to 70
GC_labels <- GC_result_hellinger_new$label_attribution
label_df <- melt(GC_labels, c("lon", "lat"), value.name = "label_attribution")
label_df$lat <- label_df$lat -70  # Adjust latitudes if necessary

# Create the plot
h <- ggplot() +
  geom_tile(data = label_df, aes(x = lon, y = lat, fill = factor(label_attribution))) +
  scale_fill_manual(values = as.vector(color_palette), na.value = NA, guide = FALSE) +  # guide = FALSE to remove legend
  ggtitle('Label attribution for GC hybrid') +
  borders("world2", colour = 'black', lwd = 0.12) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(-70, 70), expand = c(0, 0)) +  # Set y-axis limits to -70 to 70
  theme(legend.position = 'bottom') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.background = element_blank()) +
  xlab('Longitude') +
  ylab('Latitude') +
  labs(fill = 'Hellinger \nDistance') +
  theme_bw() +
  theme(legend.key.size = unit(1, 'cm'),       # Change legend key size
        legend.key.height = unit(1.4, 'cm'),   # Change legend key height
        legend.key.width = unit(0.4, 'cm'),    # Change legend key width
        legend.title = element_text(size = 16), # Change legend title font size
        legend.text = element_text(size = 12)) + # Change legend text font size
  theme(plot.title = element_text(size = 24),
        plot.subtitle = element_text(size = 20, hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) +
  easy_center_title()

h



# Generate the polychrome color palette with 26 colors
color_palette <- pals::glasbey(length(model_names))

# Adjust your labels as before
GC_labels <- GC_result_hellinger_new$label_attribution
label_df <- melt(GC_labels, c("lon", "lat"), value.name = "label_attribution")
label_df$lat <- label_df$lat - 71

# Create a spatial map using coord_sf to avoid pole artifacts
h <- ggplot() +
  geom_tile(data = label_df, aes(x = lon, y = lat, fill = factor(label_attribution))) +
  scale_fill_manual(values = as.vector(color_palette), na.value = NA, guide = FALSE) +  # Remove legend
  ggtitle('Label attribution for GC hybrid') +
  borders("world", colour = 'black', lwd = 0.12) +  # Use the "world" map instead of "world2"
  coord_sf(xlim = c(0, 359), ylim = c(-70, 70), expand = FALSE) +  # Set lat/lon limits directly
  theme(legend.position = 'bottom') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.background = element_blank()) +
  xlab('Longitude') +
  ylab('Latitude') +
  labs(fill = 'Hellinger \nDistance') +
  theme_bw() +
  theme(legend.key.size = unit(1, 'cm'),        # Change legend key size
        legend.key.height = unit(1.4, 'cm'),    # Change legend key height
        legend.key.width = unit(0.4, 'cm'),     # Change legend key width
        legend.title = element_text(size = 16), # Change legend title font size
        legend.text = element_text(size = 12)) +  # Change legend text font size
  theme(plot.title = element_text(size = 24),
        plot.subtitle = element_text(size = 20, hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) +
  easy_center_title()

h




h <- ggplot() +
  geom_tile(data = label_df, aes(x = lon, y = lat, fill = factor(label_attribution))) +
  scale_fill_manual(values = as.vector(color_palette), na.value = NA, guide = FALSE) +  # Remove legend
  ggtitle('Label attribution for GC hybrid') +
  borders("world2", colour = 'black', lwd = 0.12) +  # Use the "world2" map
  coord_quickmap(xlim = c(0, 360), ylim = c(-70, 70), expand = FALSE) +  # Use `coord_quickmap()` for better latitude/longitude handling
  theme(legend.position = 'bottom') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.background = element_blank()) +
  xlab('Longitude') +
  ylab('Latitude') +
  labs(fill = 'Hellinger \nDistance') +
  theme_bw() +
  theme(legend.key.size = unit(1, 'cm'),        # Change legend key size
        legend.key.height = unit(1.4, 'cm'),    # Change legend key height
        legend.key.width = unit(0.4, 'cm'),     # Change legend key width
        legend.title = element_text(size = 16), # Change legend title font size
        legend.text = element_text(size = 12)) +  # Change legend text font size
  theme(plot.title = element_text(size = 24),
        plot.subtitle = element_text(size = 20, hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) +
  easy_center_title()

h


h <- ggplot() +
  geom_tile(data = label_df, aes(x = lon, y = lat, fill = factor(label_attribution))) +
  scale_fill_manual(values = as.vector(color_palette), na.value = NA, guide = FALSE) +  # Remove legend
  ggtitle('Label attribution for GC hybrid') +
  borders("world2", colour = 'black', lwd = 0.12) +  # Use the "world2" map
  coord_quickmap(xlim = c(0, 360), ylim = c(-70, 70), expand = FALSE) +  # Use `coord_quickmap()` for better latitude/longitude handling
  theme(legend.position = 'bottom') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.background = element_blank()) +
  xlab('Longitude') +
  ylab('Latitude') +
  labs(fill = 'Hellinger \nDistance') +
  theme_bw() +
  theme(legend.key.size = unit(1, 'cm'),        # Change legend key size
        legend.key.height = unit(1.4, 'cm'),    # Change legend key height
        legend.key.width = unit(0.4, 'cm'),     # Change legend key width
        legend.title = element_text(size = 16), # Change legend title font size
        legend.text = element_text(size = 12)) +  # Change legend text font size
  theme(plot.title = element_text(size = 24),
        plot.subtitle = element_text(size = 20, hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) +
  easy_center_title()

h
