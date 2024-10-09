calculate_ranges <- function(variable, model_names, data_dir, year_interest, lon, lat) {
  range_var <- array(data = NA, dim = c(length(lon), length(lat), 2, length(model_names)))

  # Loop through models
  m <- 1
  for (model_name in model_names) {
    # Variable
    dir_path <- paste0(data_dir, model_name, '/', variable, '/')
    # Create the pattern
    pattern <- glob2rx(paste0(variable, "_", model_name, "*.nc"))

    # Get the filepath
    file_name <- list.files(path = dir_path, pattern = pattern)
    file_path <- paste0(dir_path, file_name)
    print(file_path)
    print(format(Sys.time(), "%Y%m%d%H%M"))

    # Check that there is only one matching file
    nc_var <<- nc_open(file_path)

    # Extract and average data for present
    yyyy <- substr(as.character(nc.get.time.series(nc_var)), 1, 4)
    iyyyy <- which(yyyy %in% year_interest)

    # Get the dimensions of the longitude and latitude variables
    lon_var <- ncvar_get(nc_var, "lon")
    lat_var <- ncvar_get(nc_var, "lat")

    # Find the indices that match the requested lon and lat ranges
    lon_indices <- which(lon_var %in% lon)
    lat_indices <- which(lat_var %in% lat)

    # Ensure the indices are contiguous ranges
    start_lon <- min(lon_indices)
    count_lon <- length(lon_indices)
    start_lat <- min(lat_indices)
    count_lat <- length(lat_indices)

    # Get the 3D data as array
    tmp_grid_var <- ncvar_get(nc_var, variable, start = c(start_lon, start_lat, min(iyyyy)),
                              count = c(count_lon, count_lat, length(iyyyy)))

    if (variable == 'pr') {
      tmp_grid_var <- log((tmp_grid_var) + 1)
    }

    range_var[ , , 1, m] <- apply(tmp_grid_var, c(1, 2), min, na.rm = TRUE)
    range_var[ , , 2, m] <- apply(tmp_grid_var, c(1, 2), max, na.rm = TRUE)

    # Close the file
    nc_close(nc_var)
    # Update counter for models
    m <- m + 1
  }

  return(range_var)
}


# Install and load necessary libraries
list_of_packages <- read.table("package_list.txt", sep="\n")$V1
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages))
  install.packages(new.packages, repos = "https://cloud.r-project.org")

library(devtools)
lapply(list_of_packages, library, character.only = TRUE)
install_github("thaos/gcoWrapR")


# Loading local functions
source_code_dir <- 'functions/' #The directory where all functions are saved.
file_paths <- list.files(source_code_dir, full.names = T)
for(path in file_paths){source(path)}


# Method 3
model_names <- read.table('model_names_long.txt')
model_names <- as.list(model_names[['V1']])

# Index of the reference
ref_index <<- 1

# Setting global variables
lon <- 0:359
lat <- -90:90
# Temporal ranges
year_interest <- 1950:2022
# data directory
data_dir <- 'data/CMIP6_merged_all/'


# List of variables
variables <- c('pr', 'tas', 'psl')

# Calculate ranges for the first variable
range_var_1 <- calculate_ranges(variables[1], model_names, data_dir, year_interest, lon, lat)
# Calculate ranges for the second variable
range_var_2 <- calculate_ranges(variables[2], model_names, data_dir, year_interest, lon, lat)
# Calculate ranges for the third variable
range_var_3 <- calculate_ranges(variables[3], model_names, data_dir, year_interest, lon, lat)

saveRDS(range_var_1, 'ranges/pr_log_range_AllModelsPar_1950-2023_90deg_3v.rds', compress = FALSE)
saveRDS(range_var_2, 'ranges/tas_range_AllModelsPar_1950-2023_90deg_3v.rds', compress = FALSE)
saveRDS(range_var_3, 'ranges/psl_range_AllModelsPar_1950-2023_90deg_3v.rds', compress = FALSE)

# Initialize the final range list for all variables
range_var_final <- list()
range_var_final[[variables[1]]] <- array(data = NA, dim = c(length(lon), length(lat), 2))
range_var_final[[variables[2]]] <- array(data = NA, dim = c(length(lon), length(lat), 2))
range_var_final[[variables[3]]] <- array(data = NA, dim = c(length(lon), length(lat), 2))

# For each grid point...
for (i in seq_along(lon)) {
  for (j in seq_along(lat)) {
    range_var_final[[variables[1]]][i, j, 1] <- min(range_var_1[i, j, 1, ])
    range_var_final[[variables[1]]][i, j, 2] <- max(range_var_1[i, j, 2, ])

    range_var_final[[variables[2]]][i, j, 1] <- min(range_var_2[i, j, 1, ])
    range_var_final[[variables[2]]][i, j, 2] <- max(range_var_2[i, j, 2, ])

    range_var_final[[variables[3]]][i, j, 1] <- min(range_var_3[i, j, 1, ])
    range_var_final[[variables[3]]][i, j, 2] <- max(range_var_3[i, j, 2, ])
  }
}

# Save the final combined ranges
saveRDS(range_var_final, 'ranges/range_var_final_allModelsPar_1950-2023_90deg_3v.rds', compress = FALSE)
