# Install and load necessary libraries
list_of_packages <- read.table("package_list.txt", sep = "\n")$V1
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages, repos = "https://cloud.r-project.org")

library(doParallel)
library(ncdf4)
library(devtools)
lapply(list_of_packages, library, character.only = TRUE)
install_github("thaos/gcoWrapR")

# Loading local functions
source_code_dir <- 'functions/'  # The directory where all functions are saved.
file_paths <- list.files(source_code_dir, full.names = TRUE)
for (path in file_paths) {
  source(path)
}

# Global variables and settings
lon <- 0:359
lat <- -70:70
year_interest <- 1950:2022
data_dir <- 'data/CMIP6_merged_all/'
variables <- c('pr', 'tas', 'psl')
model_names <- read.table('model_names_long.txt')$V1
ref_index <<- 1

# Initialize parallel backend
cl <- makeCluster(detectCores() - 1)  # Use one less core than available
registerDoParallel(cl)

# Function to calculate ranges for each variable across all models
calculate_ranges_parallel <- function(variable, model_names, data_dir, year_interest, lon, lat) {
  range_var <- array(data = NA, dim = c(length(lon), length(lat), 2, length(model_names)))  # (lon, lat, min/max, models)

  # Loop through each model (sequential)
  m <- 1
  for (model_name in model_names) {
    # File path setup
    dir_path <- paste0(data_dir, model_name, '/', variable, '/')
    pattern <- glob2rx(paste0(variable, "_", model_name, "*.nc"))
    file_name <- list.files(path = dir_path, pattern = pattern)
    file_path <- paste0(dir_path, file_name)
    print(paste("Processing model:", model_name, "with file:", file_path))
    nc_var <- nc_open(file_path)

    # Extract relevant data indices
    yyyy <- substr(as.character(nc.get.time.series(nc_var)), 1, 4)
    iyyyy <- which(yyyy %in% year_interest)
    lon_var <- ncvar_get(nc_var, "lon")
    lat_var <- ncvar_get(nc_var, "lat")
    lon_indices <- which(lon_var %in% lon)
    lat_indices <- which(lat_var %in% lat)

    # Ensure indices are contiguous ranges
    start_lon <- min(lon_indices)
    count_lon <- length(lon_indices)
    start_lat <- min(lat_indices)
    count_lat <- length(lat_indices)

    # Get 3D data as array
    tmp_grid_var <- ncvar_get(nc_var, variable, start = c(start_lon, start_lat, min(iyyyy)),
                              count = c(count_lon, count_lat, length(iyyyy)))
    if (variable == 'pr') {
      tmp_grid_var <- log(tmp_grid_var + 1)  # Apply log transformation for 'pr'
    }

    # Parallel computation of ranges over longitude and latitude points
    parallel_ranges <- foreach(i = 1:count_lon, .combine = 'rbind', .packages = 'ncdf4') %dopar% {
      local_ranges <- matrix(NA, nrow = count_lat, ncol = 2)  # Min/Max matrix for each latitude slice
      for (j in 1:count_lat) {
        local_min <- min(tmp_grid_var[i, j, ], na.rm = TRUE)
        local_max <- max(tmp_grid_var[i, j, ], na.rm = TRUE)
        local_ranges[j, ] <- c(local_min, local_max)
      }
      cbind(rep(lon_indices[i], count_lat), lat_indices, local_ranges)  # Return (lon_idx, lat_idx, min, max)
    }

    # Merge results into main range_var array
    for (k in 1:nrow(parallel_ranges)) {
      lon_idx <- parallel_ranges[k, 1]
      lat_idx <- parallel_ranges[k, 2]
      range_var[lon_idx, lat_idx, 1, m] <- parallel_ranges[k, 3]  # Min
      range_var[lon_idx, lat_idx, 2, m] <- parallel_ranges[k, 4]  # Max
    }

    # Close NetCDF file
    nc_close(nc_var)
    m <- m + 1  # Increment model counter
  }

  return(range_var)
}

# Calculate ranges for each variable using the parallelized function
range_var_1 <- calculate_ranges_parallel(variables[1], model_names, data_dir, year_interest, lon, lat)
range_var_2 <- calculate_ranges_parallel(variables[2], model_names, data_dir, year_interest, lon, lat)
range_var_3 <- calculate_ranges_parallel(variables[3], model_names, data_dir, year_interest, lon, lat)

# Save the results
saveRDS(range_var_1, 'ranges/pr_log_range_AllModelsPar_1950-2022_70deg.rds', compress = FALSE)
saveRDS(range_var_2, 'ranges/tas_range_AllModelsPar_1950-2022_70deg.rds', compress = FALSE)
saveRDS(range_var_3, 'ranges/psl_range_AllModelsPar_1950-2022_70deg.rds', compress = FALSE)

# Final global range computation by merging results
range_var_final <- list()
for (var in variables) {
  range_var_final[[var]] <- array(data = NA, dim = c(length(lon), length(lat), 2))  # Initialize final range matrix
  for (i in 1:length(lon)) {
    for (j in 1:length(lat)) {
      range_var_final[[var]][i, j, 1] <- min(get(paste0("range_var_", which(variables == var)))[i, j, 1, ], na.rm = TRUE)
      range_var_final[[var]][i, j, 2] <- max(get(paste0("range_var_", which(variables == var)))[i, j, 2, ], na.rm = TRUE)
    }
  }
}

# Save final merged ranges
saveRDS(range_var_final, 'ranges/range_var_final_allModelsPar_1950-2022_70deg.rds', compress = FALSE)

# Stop parallel backend
stopCluster(cl)
