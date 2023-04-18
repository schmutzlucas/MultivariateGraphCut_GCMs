# new color : #004DAB
h_dist_map <- array(NA, dim = dim(GC_result_hellinger$label_attribution))

for(j in seq_along(variables)){
  for(l in 0:(length(model_names))){
    islabel <- which(GC_result_hellinger$label_attribution == l)
    h_dist_map[islabel] <- h_dist_future[,,(l)][islabel]
  }
}

test_df <- melt(h_dist_map, c("lon", "lat"), value.name = "Bias")

p1 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias))+
  ggtitle(paste0('Hellinger distance of model ... ', 'Average h_dist = ', mean(h_dist_map)))+
  scale_fill_gradient(low = "white", high = "#0072B2", limits = c(0, 1))+
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
  labs(fill='Hellinger \nDistance')+
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
p1




h_dist_map <- h_dist_future[,,1]


test_df <- melt(h_dist_map, c("lon", "lat"), value.name = "Bias")

p2 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias))+
  ggtitle(paste0('Hellinger distance of model ... ', 'Average h_dist = ', mean(h_dist_map)))+
  scale_fill_gradient(low = "white", high = "#0072B2", limits = c(0, 1))+
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
  labs(fill='[K]')+
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
p2


h_dist_map <- h_dist_future[,,2]


test_df <- melt(h_dist_map, c("lon", "lat"), value.name = "Bias")

p3 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias))+
  ggtitle(paste0('Hellinger distance of model ... ', 'Average h_dist = ', mean(h_dist_map)))+
  scale_fill_gradient(low = "white", high = "#0072B2", limits = c(0, 1))+
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
  labs(fill='[K]')+
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
p3




h_dist_map <- h_dist_future[,,3]


test_df <- melt(h_dist_map, c("lon", "lat"), value.name = "Bias")

p4 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias))+
  ggtitle(paste0('Hellinger distance of model ... ', 'Average h_dist = ', mean(h_dist_map)))+
  scale_fill_gradient(low = "white", high = "#0072B2", limits = c(0, 1))+
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
  labs(fill='[K]')+
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
p4


h_dist_map <- array(NA, dim = dim(GC_result$label_attribution))

for(j in seq_along(variables)){
  for(l in 0:(length(model_names))){
    islabel <- which(GC_result$label_attribution == l)
    h_dist_map[islabel] <- h_dist_future[,,(l)][islabel]
  }
}

test_df <- melt(h_dist_map, c("lon", "lat"), value.name = "Bias")

p5 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias))+
  ggtitle(paste0('Hellinger distance of model ... ', 'Average h_dist = ', mean(h_dist_map)))+
  scale_fill_gradient(low = "white", high = "#0072B2", limits = c(0, 1))+
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
  labs(fill='[K]')+
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
p5




h_dist_map <- array(NA, dim = dim(MinBias_labels))

for(j in seq_along(variables)){
  for(l in 0:(length(model_names)-1)){
    islabel <- which(MinBias_labels == l)
    h_dist_map[islabel] <- h_dist_future[,,(l+1)][islabel]
  }
}

test_df <- melt(h_dist_map, c("lon", "lat"), value.name = "Bias")

p6 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias))+
  ggtitle(paste0('Hellinger distance of model ... ', 'Average h_dist = ', mean(h_dist_map)))+
  scale_fill_gradient(low = "white", high = "#0072B2", limits = c(0, 1))+
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
  labs(fill='[K]')+
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
p6


h_dist_map <- array(NA, dim = dim(MinBiasHellinger$label_attribution))

for(j in seq_along(variables)){
  for(l in 0:(length(model_names))){
    islabel <- which(MinBiasHellinger$label_attribution == l)
    h_dist_map[islabel] <- h_dist_future[,,(l)][islabel]
  }
}

test_df <- melt(h_dist_map, c("lon", "lat"), value.name = "Bias")

p7 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias))+
  ggtitle(paste0('Hellinger distance of model ... ', 'Average h_dist = ', mean(h_dist_map)))+
  scale_fill_gradient(low = "white", high = "#0072B2", limits = c(0, 1))+
  borders("world2", colour = 'black', lwd = 0.12) +
  scale_x_continuous(, expand = c(0, 0)) +
  scale_y_continuous(, expand = c(0,0))+
  theme(legend.position = 'bottom')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())+
  xlab('Longitude')+
  ylab('Latitude') +
  labs(fill='[K]')+
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
p7

i <- 1
h_dist_plot_models_future <- list()
for (model_name in model_names){

  h_dist_map <- h_dist_future[,,i]

  test_df <- melt(h_dist_map, c("lon", "lat"), value.name = "Bias")


  h_dist_plot_models_future[[model_name]] <- ggplot() +
    geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias))+
    labs(subtitle = 'Projection period : 1999 - 2014')+
    ggtitle(paste0(model_names[i], ': Mean Hellinger distance = ', round(mean(h_dist_map), 2)))+
    scale_fill_gradient(low = "white", high = "#015a8c", limits = c(0.1, 0.70), oob = scales::squish)+
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
  h_dist_plot_models_future[[i]]

  name <- paste0('figure/H_dist_future_',model_name)
  ggsave(paste0(name, '.pdf'), plot = h_dist_plot_models_future[[i]], width = 35, height = 25, units = "cm", dpi = 300)
  ggsave(paste0(name, '.png'), plot = h_dist_plot_models_future[[i]], width = 35, height = 25, units = "cm", dpi = 300)
  i <- i + 1
}


h_dist_map <- array(NA, dim = dim(GC_result_hellinger$label_attribution))

for(j in seq_along(variables)){
  for(l in 0:(length(model_names))){
    islabel <- which(GC_result_hellinger$label_attribution == l)
    h_dist_map[islabel] <- h_dist_future[,,(l)][islabel]
  }
}

test_df <- melt(h_dist_map, c("lon", "lat"), value.name = "Bias")

p1 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias))+
  labs(subtitle = 'Projection period : 1999 - 2014')+
  ggtitle(paste0('GraphCut Result', ': Mean Hellinger distance = ', round(mean(h_dist_map), 2)))+
  scale_fill_gradient(low = "white", high = "#015a8c", limits = c(0.1, 0.70), oob = scales::squish)+
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
p1

name <- paste0('figure/H_dist_future_GC_result')
ggsave(paste0(name, '.pdf'), plot = p1, width = 35, height = 25, units = "cm", dpi = 300)
ggsave(paste0(name, '.png'), plot = p1, width = 35, height = 25, units = "cm", dpi = 300)





h_dist_map <- array(NA, dim = dim(MinBiasHellinger$label_attribution))

for(j in seq_along(variables)){
  for(l in 0:(length(model_names))){
    islabel <- which(MinBiasHellinger$label_attribution == l)
    h_dist_map[islabel] <- h_dist_future[,,(l)][islabel]
  }
}

test_df <- melt(h_dist_map, c("lon", "lat"), value.name = "Bias")

p1 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias))+
  labs(subtitle = 'Projection period : 1999 - 2014')+
  ggtitle(paste0('GraphCut Result', ': Mean Hellinger distance = ', round(mean(h_dist_map), 2)))+
  scale_fill_gradient(low = "white", high = "#015a8c", limits = c(0.1, 0.70), oob = scales::squish)+
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
p1

name <- paste0('figure/H_dist_future_GCMinBias_result')
ggsave(paste0(name, '.pdf'), plot = p1, width = 35, height = 25, units = "cm", dpi = 300)
ggsave(paste0(name, '.png'), plot = p1, width = 35, height = 25, units = "cm", dpi = 300)