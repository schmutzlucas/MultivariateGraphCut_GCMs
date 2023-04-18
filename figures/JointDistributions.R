#############

tas_min <- min(y_breaks[18,85,,])
tas_max <- max(y_breaks[18,85,,])
tas_limit <- c(tas_min, tas_max)
pr_max <- max(x_breaks[18,85,,])
pr_limit <- c(0, pr_max)
tmp1 <- c(oTP$gc_hybrid$var[[1]])
tmp2 <- c(TPo$gc_hybrid$var[[2]])
tmp1[tmp1 < tas_min] <- NA
tmp2[tmp2 > pr_max] <- NA

title <- 'oTPo GC Hybrid 1999-2019'

df <- data.frame(x=tmp1, y=tmp2)

p <- ggplot(df, aes(x=x, y=y) ) +
  geom_bin2d(bins = 256) +
  scale_fill_gradientn(colors = viridis(256), name = 'log count', trans = 'log10', ,limits = c(NA, 31)) +
  labs(title = title) +
  xlab('Temperature [K]')+
  ylab('Precipitation [mm/day]') +
  theme_bw()+
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(0.25, 'cm'), #change legend key width
        legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=7)) + #change legend text font size
  theme(plot.title = element_text(size=16),
        axis.text=element_text(size=7),
        axis.title=element_text(size=9),)+
  scale_x_continuous(limits = c(tas_min, tas_max), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, pr_max), expand = c(0,0))+
  easy_center_title()

p





n <- 32
m <- 32

grid_data <- expand.grid(x = 1:n, y = 1:m)
grid_data$value <- pdf_matrix[18,85,,1]

grid_data$x_break <- x_breaks[18,85,,1]
grid_data$y_break <- y_breaks[18,85,,1]

ggplot(grid_data, aes(x = x_break, y = y_break, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  coord_equal() +
  theme_minimal()

pdf_matrix_2d <- matrix(pdf_matrix[18,85,,1], ncol = 32)
melted_data <- melt(pdf_matrix_2d)
colnames(melted_data) <- c("x", "y", "value")


ggplot(melted_data, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_x_continuous(breaks = 1:33, labels = x_breaks[18,85,,1]) +
  scale_y_continuous(breaks = 1:33, labels = y_breaks[18,85,,1]) +
  scale_fill_gradient(low = "white", high = "black") +
  coord_equal() +
  theme_minimal()


library(ggplot2)
library(tidyverse)
library(viridis) # for viridis color palette
library(cowplot) # for easy_center_title

# Use your own data for pdf_matrix, x_breaks, and y_breaks
pdf_matrix_2d <- matrix(pdf_matrix[18, 85, , 1], ncol = 32)
melted_data <- melt(pdf_matrix_2d)
colnames(melted_data) <- c("x", "y", "value")

# Set appropriate limits and labels for the axes
x_limits <- range(x_breaks[18, 85, , 1])
y_limits <- range(y_breaks[18, 85, , 1])
x_label <- "Temperature [K]"
y_label <- "Precipitation [mm/day]"
title <- "Your Title"

# Create the ggplot2 figure with the desired theme
ggplot(melted_data, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  theme_bw() +
  theme(legend.key.size = unit(0.5, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(0.25, 'cm'),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7),
        plot.title = element_text(size = 16),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9)) +
  easy_center_title()


########################
# Loading for models for exemples of distributions
{

  tas_min <- min(292)
  tas_max <- max(304)
  tas_limit <- c(tas_min, tas_max)
  pr_max <- max(125)
  pr_limit <- c(0, pr_max)

  lon_interest <- 18
  lat_interest <- 85

  #Color scale limits
  col_lim <- c(0,500)

  ## Model 1

  nc <- nc_open('data/CMIP6/IPSL-CM6A-LR/tas/tas_IPSL-CM6A-LR_historical_r1i1p1f1_19500101-20141230_merged_regridded_v20190614.nc')
  var <- 'tas'
  # Temporal ranges
  year_present <- 1970:2014
  yyyy <- substr(as.character(nc.get.time.series(nc)), 1, 4)
  iyyyy <- which(yyyy %in% year_present)
  tas <- ncvar_get(nc, var, start = c(lon_interest, lat_interest, min(iyyyy)), count = c(1, 1, length(iyyyy)))

  var <- 'pr'
  nc <- nc_open('data/CMIP6/IPSL-CM6A-LR/pr/pr_IPSL-CM6A-LR_historical_r1i1p1f1_19500101-20141230_merged_regridded_v20180803.nc')
  pr <- ncvar_get(nc, var, start = c(lon_interest, lat_interest, min(iyyyy)), count = c(1, 1, length(iyyyy)))
  pr <- pr * 86400



  #############



  title <- 'IPSL-CM6A-LR'

  df <- data.frame(x=tas, y=pr)

  p1 <- ggplot(df, aes(x=x, y=y) ) +
    geom_bin2d(bins = 64) +
    scale_fill_gradientn(colors = viridis(64), name = 'count', limits = col_lim, oob = scales::squish) +
    labs(title = title) +
    xlab('Temperature [K]')+
    ylab('Precipitation [mm/day]') +
    theme_bw()+
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(0.25, 'cm'), #change legend key width
          legend.title = element_text(size=9), #change legend title font size
          legend.text = element_text(size=7)) + #change legend text font size
    theme(plot.title = element_text(size=16),
          axis.text=element_text(size=7),
          axis.title=element_text(size=9),)+
    scale_x_continuous(limits = c(tas_min, tas_max), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, pr_max), expand = c(0,0))+
    easy_center_title()





  ## Model 2

  nc <- nc_open('data/CMIP6/MIROC6/tas/tas_MIROC6_historical_r1i1p1f1_19500101-20141230_merged_regridded_v20191016.nc')
  var <- 'tas'
  # Temporal ranges
  year_present <- 1970:2014
  yyyy <- substr(as.character(nc.get.time.series(nc)), 1, 4)
  iyyyy <- which(yyyy %in% year_present)
  tas <- ncvar_get(nc, var, start = c(lon_interest, lat_interest, min(iyyyy)), count = c(1, 1, length(iyyyy)))

  var <- 'pr'
  nc <- nc_open('data/CMIP6/MIROC6/pr/pr_MIROC6_historical_r1i1p1f1_19500101-20141230_merged_regridded_v20191016.nc')
  pr <- ncvar_get(nc, var, start = c(lon_interest, lat_interest, min(iyyyy)), count = c(1, 1, length(iyyyy)))
  pr <- pr * 86400



  #############


  title <- 'MIROC6'

  df <- data.frame(x=tas, y=pr)

  p2 <- ggplot(df, aes(x=x, y=y) ) +
    geom_bin2d(bins = 64) +
    scale_fill_gradientn(colors = viridis(64), name = 'count', limits = col_lim, oob = scales::squish) +
    labs(title = title) +
    xlab('Temperature [K]')+
    ylab('Precipitation [mm/day]') +
    theme_bw()+
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(0.25, 'cm'), #change legend key width
          legend.title = element_text(size=9), #change legend title font size
          legend.text = element_text(size=7)) + #change legend text font size
    theme(plot.title = element_text(size=16),
          axis.text=element_text(size=7),
          axis.title=element_text(size=9),)+
    scale_x_continuous(limits = c(tas_min, tas_max), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, pr_max), expand = c(0,0))+
    easy_center_title()



  ## Model 3

  nc <- nc_open('data/CMIP6/NorESM2-MM/tas/tas_NorESM2-MM_historical_r1i1p1f1_19500101-20141230_merged_regridded_v20191108.nc')
  var <- 'tas'
  # Temporal ranges
  year_present <- 1970:2014
  yyyy <- substr(as.character(nc.get.time.series(nc)), 1, 4)
  iyyyy <- which(yyyy %in% year_present)
  tas <- ncvar_get(nc, var, start = c(lon_interest, lat_interest, min(iyyyy)), count = c(1, 1, length(iyyyy)))

  var <- 'pr'
  nc <- nc_open('data/CMIP6/NorESM2-MM/pr/pr_NorESM2-MM_historical_r1i1p1f1_19500101-20141230_merged_regridded_v20191108.nc')
  pr <- ncvar_get(nc, var, start = c(lon_interest, lat_interest, min(iyyyy)), count = c(1, 1, length(iyyyy)))
  pr <- pr * 86400



  #############

  title <- 'NorESM2-MM'

  df <- data.frame(x=tas, y=pr)

  p3 <- ggplot(df, aes(x=x, y=y) ) +
    geom_bin2d(bins = 64) +
    scale_fill_gradientn(colors = viridis(64), name = 'count', limits = col_lim, oob = scales::squish) +
    labs(title = title) +
    xlab('Temperature [K]')+
    ylab('Precipitation [mm/day]') +
    theme_bw()+
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(0.25, 'cm'), #change legend key width
          legend.title = element_text(size=9), #change legend title font size
          legend.text = element_text(size=7)) + #change legend text font size
    theme(plot.title = element_text(size=16),
          axis.text=element_text(size=7),
          axis.title=element_text(size=9),)+
    scale_x_continuous(limits = c(tas_min, tas_max), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, pr_max), expand = c(0,0))+
    easy_center_title()




  ## Model 4

  nc <- nc_open('data/CMIP6/UKESM1-0-LL/tas/tas_UKESM1-0-LL_historical_r1i1p1f2_19500101-20141230_merged_regridded_v20190627.nc')
  var <- 'tas'
  # Temporal ranges
  year_present <- 1970:2014
  yyyy <- substr(as.character(nc.get.time.series(nc)), 1, 4)
  iyyyy <- which(yyyy %in% year_present)
  tas <- ncvar_get(nc, var, start = c(lon_interest, lat_interest, min(iyyyy)), count = c(1, 1, length(iyyyy)))

  var <- 'pr'
  nc <- nc_open('data/CMIP6/UKESM1-0-LL/pr/pr_UKESM1-0-LL_historical_r1i1p1f2_19500101-20141230_merged_regridded_v20190627.nc')
  pr <- ncvar_get(nc, var, start = c(lon_interest, lat_interest, min(iyyyy)), count = c(1, 1, length(iyyyy)))
  pr <- pr * 86400



  #############


  title <- 'UKESM1-0-LL'

  df <- data.frame(x=tas, y=pr)

  p4 <- ggplot(df, aes(x=x, y=y) ) +
    geom_bin2d(bins = 64) +
    scale_fill_gradientn(colors = viridis(64), name = 'count', limits = col_lim, oob = scales::squish) +
    labs(title = title) +
    xlab('Temperature [K]')+
    ylab('Precipitation [mm/day]') +
    theme_bw()+
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(0.25, 'cm'), #change legend key width
          legend.title = element_text(size=9), #change legend title font size
          legend.text = element_text(size=7)) + #change legend text font size
    theme(plot.title = element_text(size=16),
          axis.text=element_text(size=7),
          axis.title=element_text(size=9),)+
    scale_x_continuous(limits = c(tas_min, tas_max), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, pr_max), expand = c(0,0))+
    easy_center_title()



  # Assuming you already have p1, p2, p3, and p4 plots created
  # Modify p2, p3, and p4 to remove the legends
  p2 <- p2 + theme(legend.position = "none")
  p3 <- p3 + theme(legend.position = "none")
  p4 <- p4 + theme(legend.position = "none")

  library(cowplot)

  # Remove the legend from p1
  p1_no_legend <- p1 + theme(legend.position = "none")

  # Extract the colorbar from p1
  colorbar <- get_legend(p1)


  # Remove axis legends for the inner plots
  p1_no_legend <- p1_no_legend + theme(axis.title.x = element_blank())
  p2 <- p2 + theme(axis.title.x = element_blank())
  p2 <- p2 + theme(axis.title.y = element_blank())
  p4 <- p4 + theme(axis.title.y = element_blank())

  # ... (your existing code)

  # Combine the grid of plots with the colorbar on the right
  grid_plots <- plot_grid(p1_no_legend, p2, p3, p4, ncol = 2, nrow = 2)
  final_plot <- plot_grid(grid_plots, colorbar, ncol = 2, rel_widths = c(1, 0.1))

  # Add extra space for the title
  final_plot_with_space <- plot_grid(NULL, final_plot, ncol = 1, rel_heights = c(0.1, 1))

  # Add a general title to the final plot
  general_title <- "Comparison of Joint Distributions"
  final_plot_with_title <- ggdraw(final_plot_with_space) +
    draw_label(general_title, x = 0.5, y = 0.97, hjust = 0.5, vjust = 1, size = 24)

  final_plot_with_title

  name <- paste0('figure/distrib_plot' , lon_interest, lat_interest )
  ggsave(paste0(name, '.png'), plot = final_plot_with_title, width = 30, height = 21, units = "cm", dpi = 300)
  ggsave(paste0(name, '.pdf'), plot = final_plot_with_title, width = 30, height = 21, units = "cm", dpi = 300)
}
