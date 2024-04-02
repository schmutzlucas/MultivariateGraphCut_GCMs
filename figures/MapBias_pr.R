
# Map of the bias that works
# Creating the dataframe
# For pr
bias_tmp <- MMM$pr - reference_list$future$pr[[1]]*86400

# Ensure you are removing NA values in the computation
avg_bias <- round(mean(abs(bias_tmp), na.rm = TRUE), 2)

# Settings the limits
limit <- c(-10, 10)

# Settings the limits

v_limits <- as.numeric(format(seq(min(limit), max(limit), len=5), digits = 3))
bias_tmp[bias_tmp < min(limit)] <- min(limit)
bias_tmp[bias_tmp > max(limit)] <- max(limit)

test_df <- melt(bias_tmp, c("lon", "lat"), value.name = "Bias")


p <- ggplot()+
  geom_tile(data=test_df, aes(x=lon, y=lat, fill=Bias),)+
  ggtitle(paste0('Multi-Model Mean'))+
  labs(subtitle = paste0('Average Bias : ', avg_bias))+
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),
                       breaks = v_limits,
                       limits = limit)+
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

name <- paste0('figure/MMM_bias_pr_26models_new')
ggsave(paste0(name, '.pdf'), plot = p, width = 35, height = 25, units = "cm", dpi = 300)
ggsave(paste0(name, '.png'), plot = p, width = 35, height = 25, units = "cm", dpi = 300)


# Map of the bias that works
# Creating the dataframe
# For pr
bias_tmp <- GC_hellinger_projections$pr - reference_list$future$pr[[1]]*86400

# Ensure you are removing NA values in the computation
rm(avg_bias)
avg_bias <- round(mean(abs(bias_tmp), na.rm = TRUE), 2)

# Settings the limits
limit <- c(-10, 10)

# Settings the limits

v_limits <- as.numeric(format(seq(min(limit), max(limit), len=5), digits = 3))
bias_tmp[bias_tmp < min(limit)] <- min(limit)
bias_tmp[bias_tmp > max(limit)] <- max(limit)

test_df <- melt(bias_tmp, c("lon", "lat"), value.name = "Bias")

p <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias),)+
  ggtitle(paste0('GraphCut Hellinger'))+
  labs(subtitle = paste0('Average Bias : ', avg_bias))+
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),
                       breaks = v_limits,
                       limits = limit)+
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

name <- paste0('figure/GC_Hellinger_bias_pr_26models_new')
ggsave(paste0(name, '.pdf'), plot = p, width = 35, height = 25, units = "cm", dpi = 300)
ggsave(paste0(name, '.png'), plot = p, width = 35, height = 25, units = "cm", dpi = 300)



##### GC_hybrid


# Map of the bias that works
# Creating the dataframe
# For pr
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
# For pr

gradient_error_GC_pr <- gradient_mae(reference_list$future$pr$ERA5, GC_hellinger_projections$pr)
limit <- 2
limits <- c(0, 2)
v_limits <- as.numeric(format(seq(-limit, limit, len=6), digits = 3))
gradient_error_GC_pr[gradient_error_GC_pr < -limit] <- -limit
gradient_error_GC_pr[gradient_error_GC_pr > limit] <- limit


test_df <- melt(gradient_error_GC_pr, c("lon", "lat"), value.name = "Bias")

# Settings the limits
# limit <- max(abs(c(min(test_df$Bias), max(test_df$Bias))))
# limit <- round_any(limit, 10, floor)


p1 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias),)+
  ggtitle(paste0('GraphCut Precipitation Gradient Error'))+
  labs(subtitle = paste0('Average Error [mm/day] : ', round(mean(abs(gradient_error_GC_pr)), 2)))+
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
  labs(fill='Error \n[mm/day]')+
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



name <- paste0('figure/pr_gradient_GCMV_new')
ggsave(paste0(name, '.pdf'), plot = p1, width = 35, height = 25, units = "cm", dpi = 300)
ggsave(paste0(name, '.png'), plot = p1, width = 35, height = 25, units = "cm", dpi = 300)




# Map of the bias that works
# Creating the dataframe
# For pr

gradient_error_GC_pr <- gradient_mae(reference_list$future$pr$ERA5, MMM$pr)
limit <- 2
limits <- c(0, 2)
v_limits <- as.numeric(format(seq(-limit, limit, len=6), digits = 3))
gradient_error_GC_pr[gradient_error_GC_pr < -limit] <- -limit
gradient_error_GC_pr[gradient_error_GC_pr > limit] <- limit

test_df <- melt(gradient_error_GC_pr, c("lon", "lat"), value.name = "Bias")



p1 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias),)+
   ggtitle(paste0('MMM Precipitation Gradient Error'))+
  labs(subtitle = paste0('Average Error [mm/day] : ', round(mean(abs(gradient_error_GC_pr)), 2)))+
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
  labs(fill='Error \n[mm/day]')+
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



name <- paste0('figure/pr_gradient_MMM_new')
ggsave(paste0(name, '.pdf'), plot = p1, width = 35, height = 25, units = "cm", dpi = 300)
ggsave(paste0(name, '.png'), plot = p1, width = 35, height = 25, units = "cm", dpi = 300)



bias_tmp <- MMM$pr - reference_list$future$pr[[1]]

mean(abs(bias_tmp))

bias_tmp <- GC_hellinger_projections$pr - reference_list$future$pr[[1]]

mean(abs(bias_tmp))

