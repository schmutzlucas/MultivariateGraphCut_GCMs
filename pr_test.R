# Install and load necessary libraries
list_of_packages <- read.table("package_list.txt", sep="\n")$V1
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages))
  install.packages(new.packages, repos = "https://cran.us.r-project.org")

library(devtools)
lapply(list_of_packages, library, character.only = TRUE)
install_github("thaos/gcoWrapR")

# Import libraries
#library(c(list_of_packages))


# Loading local functions

source_code_dir <- 'functions/' #The directory where all functions are saved.
file_paths <- list.files(source_code_dir, full.names = T)
for(path in file_paths){source(path)}


data <- list()
data_matrix <- list()

# Setting global variables
lon <<- 0:359
lat <<- -90:90
# Temporal ranges
year_present <<- 1985:1999
year_future <<- 2000:2014

nc <- nc_open(paste0('data/CMIP6_merged_all/ERA5/pr/pr_ERA5_19600101-20230431_regrid.nc'))



yyyy <- substr(as.character(nc.get.time.series(nc)), 1, 4)
iyyyy <- which(yyyy %in% year_future)

tmp <- apply(
  ncvar_get(nc, 'pr', start = c(1, 1, min(iyyyy)), count = c(-1, -1, length(iyyyy))),
  1:2,
  mean
)


era5 <- tmp


nc <- nc_open(paste0('data/CMIP6_merged_all/IPSL-CM6A-LR/pr/pr_IPSL-CM6A-LR_19500101-21001230.nc'))



yyyy <- substr(as.character(nc.get.time.series(nc)), 1, 4)
iyyyy <- which(yyyy %in% year_future)

tmp <- apply(
  ncvar_get(nc, 'pr', start = c(1, 1, min(iyyyy)), count = c(-1, -1, length(iyyyy))),
  1:2,
  mean
)

ACCESS <- tmp

pr_bias <- ACCESS - era5



avg_bias <- round(mean(abs(pr_bias), na.rm = TRUE), 2)

# Settings the limits
limit <- c(-10, 10)

# Settings the limits

v_limits <- as.numeric(format(seq(min(limit), max(limit), len=5), digits = 3))
pr_bias[pr_bias < min(limit)] <- min(limit)
pr_bias[pr_bias > max(limit)] <- max(limit)

test_df <- melt(pr_bias, c("lon", "lat"), value.name = "Bias")

p <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias),)+
  ggtitle(paste0('IPSL-CM6A-LR pr Bias'))+
  labs(subtitle = paste0('Average pr Bias : ', avg_bias))+
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


name <- paste0('figure/pr_IPSL-CM6A-LR_pr_Bias')
ggsave(paste0(name, '.pdf'), plot = p, width = 35, height = 25, units = "cm", dpi = 300)
ggsave(paste0(name, '.png'), plot = p, width = 35, height = 25, units = "cm", dpi = 300)