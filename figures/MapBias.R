
# Map of the bias that works
# Creating the dataframe
# For tas
bias_tmp <- MMM$tas - reference_list$future$tas[[1]]

test_df <- melt(bias_tmp, c("lon", "lat"), value.name = "Bias")

# Settings the limits
# limit <- max(abs(c(min(test_df$Bias), max(test_df$Bias))))
# limit <- round_any(limit, 10, floor)
limit <- 8
limits <- c(-8, 8)
v_limits <- as.numeric(format(seq(-limit, limit, len=6), digits = 3))
test_df$Bias[test_df$Bias < -limit] <- -limit
test_df$Bias[test_df$Bias > limit] <- limit

p <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat, fill=Bias),)+
  ggtitle(paste0('Multi-Model Mean'))+
  labs(subtitle = paste0('Average Bias : ', round(mean(abs(bias_tmp)), 2)))+
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),
                       breaks = v_limits,
                       limits = limits)+
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
  labs(fill='Bias [K]')+
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
p
mean(abs(bias_tmp))

name <- paste0('figure/MMM_bias_tas_26models')
ggsave(paste0(name, '.pdf'), plot = p, width = 35, height = 25, units = "cm", dpi = 300)
ggsave(paste0(name, '.png'), plot = p, width = 35, height = 25, units = "cm", dpi = 300)


# Map of the bias that works
# Creating the dataframe
# For tas
bias_tmp <- GC_hellinger_projections$tas - reference_list$future$tas[[1]]

test_df <- melt(bias_tmp, c("lon", "lat"), value.name = "Bias")

# Settings the limits
# limit <- max(abs(c(min(test_df$Bias), max(test_df$Bias))))
# limit <- round_any(limit, 10, floor)

test_df$Bias[test_df$Bias < -limit] <- -limit
test_df$Bias[test_df$Bias > limit] <- limit

p <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias),)+
  ggtitle(paste0('GraphCut Hellinger'))+
  labs(subtitle = paste0('Average Bias : ', round(mean(abs(bias_tmp)), 2)))+
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),
                       breaks = v_limits,
                       limits = limits)+
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
  labs(fill='Bias [K]')+
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
p
mean(abs(bias_tmp))

name <- paste0('figure/GC_Hellinger_bias_tas_26models')
ggsave(paste0(name, '.pdf'), plot = p, width = 35, height = 25, units = "cm", dpi = 300)
ggsave(paste0(name, '.png'), plot = p, width = 35, height = 25, units = "cm", dpi = 300)



##### MMM


# Map of the bias that works
# Creating the dataframe
# For tas
bias_tmp <- mmcombi$mmm$Bias

test_df <- melt(bias_tmp, c("lon", "lat"), value.name = "Bias")

# Settings the limits
# limit <- max(abs(c(min(test_df$Bias), max(test_df$Bias))))
# limit <- round_any(limit, 10, floor)
limit <- 8
limits <- c(-8, 8)
v_limits <- as.numeric(format(seq(-limit, limit, len=6), digits = 3))
test_df$Bias[test_df$Bias < -limit] <- -limit
test_df$Bias[test_df$Bias > limit] <- limit

p1 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat, fill=Bias),)+
  ggtitle('Multi-model mean')+
  labs(subtitle = paste0('Average Bias : ', round(mean(abs(mmcombi$mmm$Bias)), 2)))+
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),
                       breaks = v_limits,
                       limits = limits)+
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
mean(abs(bias_tmp))

name <- paste0('figure/MMM_bias')
ggsave(paste0(name, '.pdf'), plot = p1, width = 35, height = 25, units = "cm", dpi = 300)
ggsave(paste0(name, '.png'), plot = p1, width = 35, height = 25, units = "cm", dpi = 300)



##### GC_hybrid


# Map of the bias that works
# Creating the dataframe
# For tas
bias_tmp <- mmcombi$gc_hybrid$Bias

test_df <- melt(bias_tmp, c("lon", "lat"), value.name = "Bias")

# Settings the limits
# limit <- max(abs(c(min(test_df$Bias), max(test_df$Bias))))
# limit <- round_any(limit, 10, floor)
limit <- 8
limits <- c(-8, 8)
v_limits <- as.numeric(format(seq(-limit, limit, len=6), digits = 3))
test_df$Bias[test_df$Bias < -limit] <- -limit
test_df$Bias[test_df$Bias > limit] <- limit

p1 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat, fill=Bias),)+
  ggtitle(paste0('GraphCut Hybrid'))+
  labs(subtitle = paste0('Average Bias : ', round(mean(abs(mmcombi$gc_hybrid$Bias)), 2)))+
  scale_fill_gradientn(colours = rev((brewer.rdbu(100))),
                       breaks = v_limits,
                       limits = limits)+
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
mean(abs(bias_tmp))


name <- paste0('figure/GC_hybrid_bias')
ggsave(paste0(name, '.pdf'), plot = p1, width = 35, height = 25, units = "cm", dpi = 300)
ggsave(paste0(name, '.png'), plot = p1, width = 35, height = 25, units = "cm", dpi = 300)




# Map of the bias that works
# Creating the dataframe
# For tas

gradient_error_GC_tas <- gradient_mae(reference_list$future$tas$ERA5, GC_hellinger_projections$tas)

test_df <- melt(gradient_error_GC_tas, c("lon", "lat"), value.name = "Bias")

# Settings the limits
# limit <- max(abs(c(min(test_df$Bias), max(test_df$Bias))))
# limit <- round_any(limit, 10, floor)
limit <- 4
limits <- c(0, 4)
v_limits <- as.numeric(format(seq(-limit, limit, len=6), digits = 3))
test_df$Bias[test_df$Bias < -limit] <- -limit
test_df$Bias[test_df$Bias > limit] <- limit

p1 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias),)+
  ggtitle(paste0('GraphCut Temperature Gradient Error'))+
  labs(subtitle = paste0('Average Error [K] : ', round(mean(abs(gradient_error_GC_tas)), 2)))+
  scale_fill_gradientn(colors = c("white", "red"),
                       breaks = v_limits,
                       limits = limits)+
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
  labs(fill='Error \n[K]')+
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



name <- paste0('figure/t_gradient_GCMV')
ggsave(paste0(name, '.pdf'), plot = p1, width = 35, height = 25, units = "cm", dpi = 300)
ggsave(paste0(name, '.png'), plot = p1, width = 35, height = 25, units = "cm", dpi = 300)




# Map of the bias that works
# Creating the dataframe
# For tas

gradient_error_GC_tas <- gradient_mae(reference_list$future$tas$ERA5, MMM$tas)

test_df <- melt(gradient_error_GC_tas, c("lon", "lat"), value.name = "Bias")

# Settings the limits
# limit <- max(abs(c(min(test_df$Bias), max(test_df$Bias))))
# limit <- round_any(limit, 10, floor)
limit <- 4
limits <- c(0, 4)
v_limits <- as.numeric(format(seq(-limit, limit, len=6), digits = 3))
test_df$Bias[test_df$Bias < -limit] <- -limit
test_df$Bias[test_df$Bias > limit] <- limit

p1 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias),)+
   ggtitle(paste0('MMM Temperature Gradient Error'))+
  labs(subtitle = paste0('Average Error [K] : ', round(mean(abs(gradient_error_GC_tas)), 2)))+
  scale_fill_gradientn(colors = c("white", "red"),
                       breaks = v_limits,
                       limits = limits)+
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
  labs(fill='Error \n[K]')+
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



name <- paste0('figure/t_gradient_MMM')
ggsave(paste0(name, '.pdf'), plot = p1, width = 35, height = 25, units = "cm", dpi = 300)
ggsave(paste0(name, '.png'), plot = p1, width = 35, height = 25, units = "cm", dpi = 300)



bias_tmp <- MMM$pr - reference_list$future$pr[[1]]

mean(abs(bias_tmp))

bias_tmp <- GC_hellinger_projections$pr - reference_list$future$pr[[1]]

mean(abs(bias_tmp))

