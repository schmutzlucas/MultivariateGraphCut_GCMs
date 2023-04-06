library(stringr)

# Get model names in Data
dir_path <- "data/CMIP6"
model_names <- dir(dir_path)

variables <- c('tas', 'pr', 'tasmax')
exp <- 'historical'


nc <- nc_open('data/IPSL-CM6A-LR/tas/tas_IPSL-CM6A-LR_historical_r1i1p1f1_19500101-20141230_merged_regridded_v20190614.nc')
var <- 'tas'
# Temporal ranges
year_present <- 1970:1990
year_future <- 1999:2019
yyyy <- substr(as.character(nc.get.time.series(nc)), 1, 4)
iyyyy <- which(yyyy %in% year_present)
tas <- ncvar_get(nc, var, start = c(1, 1, min(iyyyy)), count = c(-1, -1, length(iyyyy)))

var <- 'pr'
nc <- nc_open('data/IPSL-CM6A-LR/pr/pr_IPSL-CM6A-LR_historical_r1i1p1f1_19500101-20141230_merged_regridded_v20180803.nc')
pr <- ncvar_get(nc, var, start = c(1, 1, min(iyyyy)), count = c(-1, -1, length(iyyyy)))

var <- 'tasmax'
nc <- nc_open('data/IPSL-CM6A-LR/tasmax/tasmax_IPSL-CM6A-LR_historical_r1i1p1f1_19500101-20141230_merged_regridded_v20190614.nc')
tas_max <- ncvar_get(nc, var, start = c(1, 1, min(iyyyy)), count = c(-1, -1, length(iyyyy)))



# Parallel backend
n_cores <- detectCores()
cl <- makeCluster(n_cores)
registerDoParallel(6)

range_pr <- c(0, 5e-3)
range_tas <- c(170, 330)
nbins <- 20
pdf_values <- array(dim = c(dim(tas)[1], dim(tas)[2], nbins^2))

# Create a function for the inner loop to simplify the code
inner_loop <- function(i, j, tas, pr, nbins, range_tas, range_pr) {
  dump <- kde2d(tas[i, j,], pr[i, j,],
                n = nbins, h = 1, lims = c(range_tas, range_pr))
  return(c(dump$z))
}

# Parallelize the outer loop
results <- foreach(i = 1:dim(tas)[1], .combine = 'rbind') %dopar% {
  inner_results <- foreach(j = 1:dim(tas)[2], .combine = 'rbind') %do% {
    inner_loop(i, j, tas, pr, nbins, range_tas, range_pr)
  }
  inner_results
}


range_pr <- c(0, 5e-3)
range_tas <- c(170, 330)
nbins <- 50
pdf_values <- array(dim = c(dim(tas)[1], dim(tas)[2], nbins^2))
for (i in 1:360) {
  for (j in 1:181) {
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

