# Generate the polychrome color palette with 26 colors
color_palette <- pals::glasbey(length(model_names))

GC_labels <- GC_result_hellinger_new$label_attribution

label_df <- melt(GC_labels, c("lon", "lat"), value.name = "label_attribution")
label_df$lat <- label_df$lat - 71

h <- ggplot() +
  geom_tile(data = label_df, aes(x = lon, y = lat, fill = factor(label_attribution))) +
  scale_fill_manual(values = as.vector(color_palette), na.value = NA, guide = FALSE) +  # guide = FALSE to remove legend
  ggtitle('Label attribution for GC hybrid') +
  borders("world2", colour = 'black', lwd = 0.12) +
  scale_x_continuous(, expand = c(0, 0)) +
  scale_y_continuous(, expand = c(0,0))+
  theme(legend.position = 'bottom')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())+
  xlab('Longitude')+
  ylab('Latitude') +
  labs(fill='Hellinger \nDistance')+
  theme_bw()+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=12))+ #change legend text font size
  theme(plot.title = element_text(size=24),
        plot.subtitle = element_text(size = 20,hjust=0.5),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),)+
  easy_center_title()

h


# Generate the polychrome color palette with 26 colors
color_palette <- pals::glasbey(26)

GC_labels <- GC_result_hellinger_new$label_attribution

label_df <- melt(GC_labels, c("lon", "lat"), value.name = "label_attribution")
label_df$lat <- label_df$lat - 91

h <- ggplot() +
  geom_tile(data = label_df, aes(x = lon, y = lat, fill = factor(label_attribution))) +
  scale_fill_manual(values = as.vector(color_palette), na.value = NA, guide = FALSE) +  # guide = FALSE to remove legend
  ggtitle('Label attribution for GC hybrid') +
  borders("world2", colour = 'black', lwd = 0.12) +
  scale_x_continuous(, expand = c(0, 0)) +
  scale_y_continuous(, expand = c(0,0))+
  theme(legend.position = 'bottom')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())+
  xlab('Longitude')+
  ylab('Latitude') +
  labs(fill='Hellinger \nDistance')+
  theme_bw()+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=12))+ #change legend text font size
  theme(plot.title = element_text(size=24),
        plot.subtitle = element_text(size = 20,hjust=0.5),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),)+
  easy_center_title()

h

name <- paste0('figure/labelling_2')
ggsave(paste0(name, '.pdf'), plot = h, width = 35, height = 25, units = "cm", dpi = 300)
ggsave(paste0(name, '.png'), plot = h, width = 35, height = 25, units = "cm", dpi = 300)



# Define the color mapping for label_attribution
# Define the color mapping for label_attribution
color_mapping <- c("1" = "#4f6980", "2" = "#bfbb60", "3" = "#638b66", "4" = "#b34f6a")



# Add a new column in label_df with the corresponding colors
label_df$color <- color_mapping[as.character(label_df$label_attribution)]

h <- ggplot() +
  geom_tile(data=label_df, aes(x=lon, y=lat, fill= color),)+
  scale_fill_identity(na.value = NA,
                      guide = guide_legend(override.aes = list(fill = color_mapping),
                                           title = "Model",
                                           labels = model_names)) +
  ggtitle('Label attribution for GC hybrid')+
  borders("world2", colour = 'black', lwd = 0.12) +
  scale_x_continuous(, expand = c(0, 0)) +
  scale_y_continuous(, expand = c(0,0))+
  theme(legend.position = 'bottom')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())+
  xlab('Longitude')+
  ylab('Latitude') +
  labs(fill='[mm/day]')+
  theme_bw()+
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.7, 'cm'), #change legend key height
        legend.key.width = unit(0.2, 'cm'), #change legend key width
        legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=7))+ #change legend text font size
  theme(plot.title = element_text(size=16),
        axis.text=element_text(size=7),
        axis.title=element_text(size=9),)+
  easy_center_title()
h
name <- paste0('figure/labelling' , lon, lat )
ggsave(paste0(name, '.png'), plot = h, width = 30, height = 21, units = "cm", dpi = 300)
ggsave(paste0(name, '.pdf'), plot = h, width = 30, height = 21, units = "cm", dpi = 300)




GC_labels <- GC_result$label_attribution

label_df <-melt(GC_labels, c("lon", "lat"),
                value.name = "label_attribution")
label_df$lat <- label_df$lat -91


p <- ggplot() +
  geom_tile(data=label_df, aes(x=lon, y=lat, fill= factor(label_attribution)),)+
  scale_fill_hue(na.value=NA,
                 labels = model_names,
                 name = 'Model')+
  ggtitle('Label attribution for GC hybrid')+
  # guide = guide_colourbar(
  #   barwidth = 0.4,
  #   ticks.colour = 'black',
  #   ticks.linewidth = 1,
  #   frame.colour = 'black',
  #   draw.ulim = TRUE,
  #   draw.llim = TRUE),)

  borders("world2", colour = 'black', lwd = 0.12) +
  scale_x_continuous(, expand = c(0, 0)) +
  scale_y_continuous(, expand = c(0,0))+
  theme(legend.position = 'bottom')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())+
  xlab('Longitude')+
  ylab('Latitude') +
  labs(fill='[mm/day]')+
  theme_bw()+
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.7, 'cm'), #change legend key height
        legend.key.width = unit(0.2, 'cm'), #change legend key width
        legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=7))+ #change legend text font size
  theme(plot.title = element_text(size=16),
        axis.text=element_text(size=7),
        axis.title=element_text(size=9),)+
  easy_center_title()


label_df <-melt(MinBias_labels, c("lon", "lat"),
                value.name = "label_attribution")
label_df$lat <- label_df$lat -91
q <- ggplot() +
  geom_tile(data=label_df, aes(x=lon, y=lat, fill= factor(label_attribution)),)+
  scale_fill_hue(na.value=NA,
                 labels = model_names,
                 name = 'Model')+
  ggtitle('Label attribution for GC hybrid')+
  # guide = guide_colourbar(
  #   barwidth = 0.4,
  #   ticks.colour = 'black',
  #   ticks.linewidth = 1,
  #   frame.colour = 'black',
  #   draw.ulim = TRUE,
  #   draw.llim = TRUE),)

  borders("world2", colour = 'black', lwd = 0.12) +
  scale_x_continuous(, expand = c(0, 0)) +
  scale_y_continuous(, expand = c(0,0))+
  theme(legend.position = 'bottom')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())+
  xlab('Longitude')+
  ylab('Latitude') +
  labs(fill='[mm/day]')+
  theme_bw()+
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.7, 'cm'), #change legend key height
        legend.key.width = unit(0.2, 'cm'), #change legend key width
        legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=7))+ #change legend text font size
  theme(plot.title = element_text(size=16),
        axis.text=element_text(size=7),
        axis.title=element_text(size=9),)+
  easy_center_title()
q



label_df <-melt(MinBiasHellinger$label_attribution, c("lon", "lat"),
                value.name = "label_attribution")
label_df$lat <- label_df$lat -91
q <- ggplot() +
  geom_tile(data=label_df, aes(x=lon, y=lat, fill= factor(label_attribution)),)+
  scale_fill_hue(na.value=NA,
                 labels = model_names,
                 name = 'Model')+
  ggtitle('Label attribution for GC hybrid')+
  # guide = guide_colourbar(
  #   barwidth = 0.4,
  #   ticks.colour = 'black',
  #   ticks.linewidth = 1,
  #   frame.colour = 'black',
  #   draw.ulim = TRUE,
  #   draw.llim = TRUE),)

  borders("world2", colour = 'black', lwd = 0.12) +
  scale_x_continuous(, expand = c(0, 0)) +
  scale_y_continuous(, expand = c(0,0))+
  theme(legend.position = 'bottom')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())+
  xlab('Longitude')+
  ylab('Latitude') +
  labs(fill='[mm/day]')+
  theme_bw()+
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.7, 'cm'), #change legend key height
        legend.key.width = unit(0.2, 'cm'), #change legend key width
        legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=7))+ #change legend text font size
  theme(plot.title = element_text(size=16),
        axis.text=element_text(size=7),
        axis.title=element_text(size=9),)+
  easy_center_title()
q
