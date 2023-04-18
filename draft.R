library(stringr)

# Get model names in Data
dir_path <- "data/CMIP6"
model_names <- dir(dir_path)

variables <- c('tas', 'pr', 'tasmax')
exp <- 'historical'


nc <- nc_open('data/CMIP6/IPSL-CM6A-LR/tas/tas_IPSL-CM6A-LR_historical_r1i1p1f1_19500101-20141230_merged_regridded_v20190614.nc')
var <- 'tas'
# Temporal ranges
year_present <- 1970:1990
year_future <- 1999:2019
yyyy <- substr(as.character(nc.get.time.series(nc)), 1, 4)
iyyyy <- which(yyyy %in% year_present)
tas <- ncvar_get(nc, var, start = c(1, 1, min(iyyyy)), count = c(1, 1, length(iyyyy)))

var <- 'pr'
nc <- nc_open('data/CMIP6/IPSL-CM6A-LR/pr/pr_IPSL-CM6A-LR_historical_r1i1p1f1_19500101-20141230_merged_regridded_v20180803.nc')
pr <- ncvar_get(nc, var, start = c(1, 1, min(iyyyy)), count = c(-1, -1, length(iyyyy)))

var <- 'tasmax'
nc <- nc_open('data/CMIP6/IPSL-CM6A-LR/tasmax/tasmax_IPSL-CM6A-LR_historical_r1i1p1f1_19500101-20141230_merged_regridded_v20190614.nc')
tas_max <- ncvar_get(nc, var, start = c(1, 1, min(iyyyy)), count = c(-1, -1, length(iyyyy)))





range_pr <- c(0, 5e-3)
range_tas <- c(170, 330)
nbins <- 50
pdf_values <- array(dim = c(dim(tas)[1], dim(tas)[2], nbins^2))
for (i in 1:length(dim(tas[1]))) {
  for (j in 1:length(dim(tas[2]))) {
    print(c(i,j))
    range_tas <- range(tas[i,j,])
    range_pr <- range(pr[i,j,])
    # times_series <- matrix(data = NA, nrow = dim(tas)[3], ncol = 2)
    # times_series[,1] <- tas[i,j,]
    # times_series[,2] <- pr[i,j,]

    dump <- kde2d(tas[i,j,], pr[i,j,],
                    n = nbins, lims = c(range_tas, range_pr))
    pdf_values[i,j,] <- c(dump$z)


  }
}

nbins <- 50
range_tas <- range(tas[45,45,])
range_pr <- range((log(pr[45,45,]*86400+1)))
dump <- kde2d(tas[45,45,], (log(pr[45,45,]*86400+1)),
                n = nbins, lims = c(range_tas, range_pr))

contour(dump)

tp_dens_model3_long <- reshape2::melt(dump$z)

# Create a ggplot graph of the density estimates
ggplot(tp_dens_model3_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_viridis() +
  labs(x = "Temperature", y = "Precipitation", fill = "Density") +
  ggtitle("Density Estimation for Climate Model 3") +
  theme_bw()


times_series <- matrix(data = NA, nrow = dim(tas)[3], ncol = 2)

times_series[,1] <- tas[100,100,]
times_series[,2] <- pr[100,100,]
nbins <- 50

joint_dist <- kde2d(times_series[, 1], times_series[, 2],
                    n = nbins, lims = c(range_tas, range_pr))

dir_path <- paste0('data/CMIP6/')
var <- 'pr'
model_name <- 'IPSL-CM6A-LR'
experiment <- 'historical'
# Create the pattern
pattern <- glob2rx(paste0(var, "_", model_name, "_", experiment, "*.nc"))

filepath <- list.files(path = dir_path,
                       pattern = pattern)



# Set the variables
var <- "pr"
model_name <- "IPSL-CM6A-LR"
experiment <- "historical"

# Create the pattern using glob2rx()
pattern <- glob2rx(paste0(var, "_", model_name, "_", experiment, "*.nc"))

# Test the pattern
test_vec <- c("pr_IPSL-CM6A-LR_historical_r1i1p1f1_19500101-20141230_merged_regridded_v20180803.nc",
              "pr_IPSL-CM6A-LR_ssp585_r1i1p1f1_20150101-21001230_merged_regridded_v20190903.nc",
              "tasmin_IPSL-CM6A-LR_historical_r1i1p1f1_19500101-20141230_merged_regridded_v20180803.nc",
              "tasmax_NorESM2-MM_historical_r1i1p1f1_19500101-20141230_merged_regridded_v20191108.nc")
grep(pattern, test_vec)

model_name <- 'MIROC6'
var <- 'tas'
period <- 'historical'

dir_path <- paste0('data/CMIP6/', model_name, '/', var, '/')
# Create the pattern
pattern <- glob2rx(paste0(var, "_", model_name, "_", period, "*.nc"))

# Get the filepath
file_path <- list.files(path = dir_path,
                       pattern = pattern)



# Generate some sample data
x <- rnorm(100)

# Compute the KDE for the range -3 to 3
dens <- density(x, from = -3, to = 3)

# Plot the density estimate
plot(dens)





nc2 <- nc_open('data/CMIP6/IPSL-CM6A-LR/tas/tas_IPSL-CM6A-LR_historical_r1i1p1f1_19500101-20141230_merged_regridded_v20190614.nc')
var <- 'tas'
# Temporal ranges
year_present <- 1970:1990
year_future <- 1999:2019
yyyy <- substr(as.character(nc.get.time.series(nc2)), 1, 4)
iyyyy <- which(yyyy %in% year_present)
tmp2 <- ncvar_get(nc2, var, start = c(1, 1, min(iyyyy)), count = c(-1, -1, length(iyyyy)))

var <- 'pr'
nc1 <- nc_open('data/CMIP6/IPSL-CM6A-LR/pr/pr_IPSL-CM6A-LR_historical_r1i1p1f1_19500101-20141230_merged_regridded_v20180803.nc')
tmp1 <- ncvar_get(nc1, var, start = c(1, 1, min(iyyyy)), count = c(-1, -1, length(iyyyy)))

vec1 <- tmp1[123, 45, ]
vec2 <- tmp2[123, 45, ]

breaks1 <- seq(min(tmp1[123, 45, ]), max(tmp1[123, 45, ]), length.out = 50) # adjust the number of breaks as needed
breaks2 <- seq(min(tmp2[123, 45, ]), max(tmp2[123, 45, ]), length.out = 50)



library(squash)
terst <- hist2(vec1, vec2, xbreaks = breaks1, ybreaks = breaks2, plot = FALSE)





  title <- 'IPSL-CM6A-LR'

  df <- data.frame(x=tas, y=pr+1)

  p1 <- ggplot(df, aes(x=x, y=y) ) +
    geom_bin2d(bins = 64) +
    scale_fill_gradientn(colors = viridis(64), name = 'count', limits = col_lim) +
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
    scale_y_log2(range = c(0, pr_max))+
    easy_center_title()

p1