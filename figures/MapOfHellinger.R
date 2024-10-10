
h_dist_map <- array(NA, dim = dim(GC_result_hellinger_new$label_attribution))

for(j in seq_along(variables)){
  for(l in 0:(length(model_names))){
    islabel <- which(GC_result_hellinger$label_attribution == l)
    h_dist_map[islabel] <- GC_result_hellinger[['h_dist']][,,(l)][islabel]
  }
}

test_df <- melt(h_dist_map, c("lon", "lat"), value.name = "Bias")

p1 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias))+
  ggtitle(paste0('Hellinger distance of model ... ', 'Average h_dist = ', mean(h_dist_map)))+
    scale_fill_gradient(low = "white", high = "blue") +
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
p1




h_dist_map <- GC_result_hellinger$h_dist[,,1]


test_df <- melt(h_dist_map, c("lon", "lat"), value.name = "Bias")

p2 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias))+
  ggtitle(paste0('Hellinger distance of model ... ', 'Average h_dist = ', mean(h_dist_map)))+
    scale_fill_gradient(low = "white", high = "blue") +
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


h_dist_map <- GC_result_hellinger$h_dist[,,2]


test_df <- melt(h_dist_map, c("lon", "lat"), value.name = "Bias")

p3 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias))+
  ggtitle(paste0('Hellinger distance of model ... ', 'Average h_dist = ', mean(h_dist_map)))+
    scale_fill_gradient(low = "white", high = "blue") +
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




h_dist_map <- GC_result_hellinger$h_dist[,,3]


test_df <- melt(h_dist_map, c("lon", "lat"), value.name = "Bias")

p4 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias))+
  ggtitle(paste0('Hellinger distance of model ... ', 'Average h_dist = ', mean(h_dist_map)))+
    scale_fill_gradient(low = "white", high = "blue") +
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
    h_dist_map[islabel] <- GC_result_hellinger[['h_dist']][,,(l)][islabel]
  }
}

test_df <- melt(h_dist_map, c("lon", "lat"), value.name = "Bias")

p5 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias))+
  ggtitle(paste0('Hellinger distance of model ... ', 'Average h_dist = ', mean(h_dist_map)))+
    scale_fill_gradient(low = "white", high = "blue") +
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
    h_dist_map[islabel] <- GC_result_hellinger[['h_dist']][,,(l+1)][islabel]
  }
}

test_df <- melt(h_dist_map, c("lon", "lat"), value.name = "Bias")

p6 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias))+
  ggtitle(paste0('Hellinger distance of model ... ', 'Average h_dist = ', mean(h_dist_map)))+
    scale_fill_gradient(low = "white", high = "blue") +
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
    h_dist_map[islabel] <- MinBiasHellinger[['h_dist']][,,(l)][islabel]
  }
}

test_df <- melt(h_dist_map, c("lon", "lat"), value.name = "Bias")

p7 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias))+
  ggtitle(paste0('Hellinger distance of model ... ', 'Average h_dist = ', mean(h_dist_map)))+
    scale_fill_gradient(low = "white", high = "blue") +
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
p7




h_dist_map <- array(NA, dim = dim(MinBias_labels))

for(j in seq_along(variables)){
  for(l in 0:(length(model_names)-1)){
    islabel <- which(MinBias_labels == l)
    h_dist_map[islabel] <- GC_result_hellinger_new[['h_dist']][,,(l+1)][islabel]
  }
}

test_df <- melt(h_dist_map, c("lon", "lat"), value.name = "Bias")

p6 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias))+
  ggtitle(paste0('Hellinger distance of model ... ', 'Average h_dist = ', mean(h_dist_map)))+
    scale_fill_gradient(low = "white", high = "blue") +
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
