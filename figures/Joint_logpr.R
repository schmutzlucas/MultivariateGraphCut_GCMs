########################
# Loading for models for exemples of distributions
{


  tas_min <- min(292)
  tas_max <- max(305)
  tas_limit <- c(tas_min, tas_max)
  pr_max <- max(log2(180))
  pr_limit <- c(0, pr_max)

  # #Vienna
  # lon <- 16
  # lat <- 138

  # # Lausanne
  # lon <- 6
  # lat <- 90 + 46

  # # Idea 1
  # lon <- 18
  # lat <- 90 -12
  #
  #   # Idea 1
  # lon <- 95
  # lat <- 90 +27

    #   # Idea 1
   lon <- 114
   lat <- 94





  #Color scale limits
  col_lim <- c(0,80)

  ## Model 1

  nc <- nc_open('data/CMIP6/FGOALS-f3-L/tas/tas_FGOALS-f3-L_historical_r1i1p1f1_19500101-20141230_merged_regridded_v20191019.nc')
  var <- 'tas'
  # Temporal ranges
  year_present <- 1970:2014
  yyyy <- substr(as.character(nc.get.time.series(nc)), 1, 4)
  iyyyy <- which(yyyy %in% year_present)
  tas <- ncvar_get(nc, var, start = c(lon, lat, min(iyyyy)), count = c(1, 1, length(iyyyy)))

  var <- 'pr'
  nc <- nc_open('data/CMIP6/FGOALS-f3-L/pr/pr_FGOALS-f3-L_historical_r1i1p1f1_19500101-20141230_merged_regridded_v20191019.nc')
  pr <- ncvar_get(nc, var, start = c(lon, lat, min(iyyyy)), count = c(1, 1, length(iyyyy)))
  pr <- log2(pr * 86400 + 1)



  #############



  title <- 'FGOALS-f3-L'

  df <- data.frame(x=tas, y=pr)

  p2 <- ggplot(df, aes(x=x, y=y))+
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
  tas <- ncvar_get(nc, var, start = c(lon, lat, min(iyyyy)), count = c(1, 1, length(iyyyy)))

  var <- 'pr'
  nc <- nc_open('data/CMIP6/MIROC6/pr/pr_MIROC6_historical_r1i1p1f1_19500101-20141230_merged_regridded_v20191016.nc')
  pr <- ncvar_get(nc, var, start = c(lon, lat, min(iyyyy)), count = c(1, 1, length(iyyyy)))
  pr <- log2(pr * 86400 + 1)



  #############


  title <- 'MIROC6'

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



  ## Model 3

  nc <- nc_open('data/CMIP6/MPI-ESM1-2-HR/tas/tas_MPI-ESM1-2-HR_historical_r1i1p1f1_19500101-20141230_merged_regridded_v20190710.nc')
  var <- 'tas'
  # Temporal ranges
  year_present <- 1970:2014
  yyyy <- substr(as.character(nc.get.time.series(nc)), 1, 4)
  iyyyy <- which(yyyy %in% year_present)
  tas <- ncvar_get(nc, var, start = c(lon, lat, min(iyyyy)), count = c(1, 1, length(iyyyy)))

  var <- 'pr'
  nc <- nc_open('data/CMIP6/MPI-ESM1-2-HR/pr/pr_MPI-ESM1-2-HR_historical_r1i1p1f1_19500101-20141230_merged_regridded_v20190710.nc')
  pr <- ncvar_get(nc, var, start = c(lon, lat, min(iyyyy)), count = c(1, 1, length(iyyyy)))
  pr <- log2(pr * 86400 + 1)



  #############

  title <- 'MPI-ESM1-2-HR'

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

  nc <- nc_open('data/CMIP6/INM-CM5-0/tas/tas_INM-CM5-0_historical_r1i1p1f1_19500101-20141230_merged_regridded_v20190610.nc')
  var <- 'tas'
  # Temporal ranges
  year_present <- 1970:2014
  yyyy <- substr(as.character(nc.get.time.series(nc)), 1, 4)
  iyyyy <- which(yyyy %in% year_present)
  tas <- ncvar_get(nc, var, start = c(lon, lat, min(iyyyy)), count = c(1, 1, length(iyyyy)))

  var <- 'pr'
  nc <- nc_open('data/CMIP6/INM-CM5-0/pr/pr_INM-CM5-0_historical_r1i1p1f1_19500101-20141230_merged_regridded_v20190610.nc')
  pr <- ncvar_get(nc, var, start = c(lon, lat, min(iyyyy)), count = c(1, 1, length(iyyyy)))
  pr <- log2(pr * 86400 + 1)



  #############


  title <- 'INM-CM5-0'

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


  library(scales)

  # ... (your existing code)



  # Remove the legend from p1
  p1_no_legend <- p1 + theme(legend.position = "none")

  # Extract the colorbar from p1
  colorbar <- get_legend(p1)


  # Remove axis legends for the inner plots
  p1_no_legend <- p1_no_legend + theme(axis.title.x = element_blank())
  p2 <- p2 + theme(axis.title.x = element_blank())
  p2 <- p2 + theme(axis.title.y = element_blank())
  p4 <- p4 + theme(axis.title.y = element_blank())



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
  name <- paste0('figure/distrib_plot_log' , lon, lat )
  ggsave(paste0(name, '.png'), plot = final_plot_with_title, width = 30, height = 21, units = "cm", dpi = 300)
  ggsave(paste0(name, '.pdf'), plot = final_plot_with_title, width = 30, height = 21, units = "cm", dpi = 300)
}

final_plot_with_title
