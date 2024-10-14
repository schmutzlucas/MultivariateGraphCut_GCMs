# -----------------------------------------------
# Script to Calculate and Merge Variable Ranges
# across Multiple Climate Models Using Parallel
# Processing in R
# -----------------------------------------------

# This script performs the following operations:
# 1. Installs and loads the necessary R packages.
# 2. Loads custom functions from local directories.
# 3. Defines global variables such as longitude, latitude, and models.
# 4. Sets up a parallel backend to speed up the processing.
# 5. Defines a function to compute the minimum and maximum range
#    values for each variable (pr, tas, psl) across a list of
#    climate models.
# 6. Uses parallel processing (`future_lapply`) to compute the ranges.
# 7. Merges the computed ranges into a final global range for each variable.
# 8. Saves the results as RDS files.
# 9. Displays the total execution time.
# -----------------------------------------------

# 1. Install and load necessary libraries
# Read the list of required packages from a text file and install if not already present
list_of_packages <- read.table("package_list.txt", sep = "\n")$V1
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages, repos = "https://cloud.r-project.org")

# Load necessary libraries
library(future)
library(future.apply)
library(ncdf4)
library(devtools)

# Load other packages mentioned in the list
lapply(list_of_packages, library, character.only = TRUE)

# Install and load a custom package for Graph Cut Optimization
install_github("thaos/gcoWrapR")

# 2. Load custom functions from a local directory
# Define the directory containing the custom functions
source_code_dir <- 'functions/'
file_paths <- list.files(source_code_dir, full.names = TRUE)
# Load all R scripts in the specified directory
for (path in file_paths) source(path)

# Start timing the script execution
start_time <- Sys.time()
cat("Script started at: ", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")

# 3. Global variables definition
# Define longitude, latitude, years of interest, data directory, variables, and model names
lon <- 0:359
lat <- -70:70
year_interest <- 1960:2022
data_dir <- 'data/CMIP6_merged_all/'
variables <- c('pr', 'tas', 'psl')
model_names <- read.table('model_names_pr_tas_psl.txt')$V1

# 4. Set up parallel backend using `future` package
# This enables parallel processing with 8 workers for faster computation
plan(multisession, workers = 4)

# 5. Function to calculate the minimum and maximum ranges for each variable
# across all models, for each pixel in the data grid.
# Inputs:
# - variable: Climate variable (e.g., pr, tas, psl)
# - model_names: List of model names
# - data_dir: Directory containing the model data
# - year_interest: Years of interest for analysis
# - lon, lat: Longitude and Latitude grids
calculate_ranges <- function(variable, model_names, data_dir, year_interest, lon, lat) {
  # Create an empty array to store the range values
  range_var <- array(NA, dim = c(length(lon), length(lat), 2, length(model_names)))

  # Loop through each model to calculate the range for the specified variable
  for (m in seq_along(model_names)) {
    model_name <- model_names[m]
    # Construct the file path and open the NetCDF file
    file_path <- paste0(data_dir, model_name, '/', variable, '/', list.files(path = paste0(data_dir, model_name, '/', variable, '/'), pattern = glob2rx(paste0(variable, "_", model_name, "*.nc")))[1])
    nc_var <- nc_open(file_path)

    # Extract time series and identify indices for the years of interest
    yyyy <- substr(as.character(nc.get.time.series(nc_var)), 1, 4)
    iyyyy <- which(yyyy %in% year_interest)

    # Get the longitude and latitude indices for the region of interest
    lon_var <- ncvar_get(nc_var, "lon")
    lat_var <- ncvar_get(nc_var, "lat")
    lon_indices <- which(lon_var %in% lon)
    lat_indices <- which(lat_var %in% lat)
    start_lon <- min(lon_indices)
    start_lat <- min(lat_indices)

    # Extract the variable data for the region and years of interest
    tmp_grid_var <- ncvar_get(nc_var, variable, start = c(start_lon, start_lat, min(iyyyy)), count = c(length(lon_indices), length(lat_indices), length(iyyyy)))

    # Transform the variable data if necessary (e.g., log transformation for precipitation)
    if (variable == 'pr') tmp_grid_var <- log(tmp_grid_var + 1)

    # Calculate the min and max ranges for each grid cell
    range_var[ , , 1, m] <- apply(tmp_grid_var, c(1, 2), min, na.rm = TRUE)
    range_var[ , , 2, m] <- apply(tmp_grid_var, c(1, 2), max, na.rm = TRUE)

    # Close the NetCDF file
    nc_close(nc_var)
  }

  return(range_var)
}

# 6. Calculate ranges for each variable asynchronously using `future_lapply`
range_results <- future_lapply(variables, function(var) calculate_ranges(var, model_names, data_dir, year_interest, lon, lat))

# 7. Assign the results to respective variables for clarity
range_var_1 <- range_results[[1]]
range_var_2 <- range_results[[2]]
range_var_3 <- range_results[[3]]

# 8. Save each variable's range separately as RDS files
saveRDS(range_var_1, 'ranges/pr_log_range_AllModelsPar_1950-2023_90deg_3v.rds', compress = FALSE)
saveRDS(range_var_2, 'ranges/tas_range_AllModelsPar_1950-2023_90deg_3v.rds', compress = FALSE)
saveRDS(range_var_3, 'ranges/psl_range_AllModelsPar_1950-2023_90deg_3v.rds', compress = FALSE)

# 9. Final global range computation by merging results for each variable
range_var_final <- list()
for (var in variables) range_var_final[[var]] <- array(NA, dim = c(length(lon), length(lat), 2))

# Loop through each grid cell to merge the results across all models
for (i in seq_along(lon)) {
  for (j in seq_along(lat)) {
    range_var_final[[variables[1]]][i, j, 1] <- min(range_var_1[i, j, 1, ], na.rm = TRUE)
    range_var_final[[variables[1]]][i, j, 2] <- max(range_var_1[i, j, 2, ], na.rm = TRUE)
    range_var_final[[variables[2]]][i, j, 1] <- min(range_var_2[i, j, 1, ], na.rm = TRUE)
    range_var_final[[variables[2]]][i, j, 2] <- max(range_var_2[i, j, 2, ], na.rm = TRUE)
    range_var_final[[variables[3]]][i, j, 1] <- min(range_var_3[i, j, 1, ], na.rm = TRUE)
    range_var_final[[variables[3]]][i, j, 2] <- max(range_var_3[i, j, 2, ], na.rm = TRUE)
  }
}

# 10. Save the final merged ranges
saveRDS(range_var_final, 'ranges/range_var_final_allModelsPar_1950-2023_90deg_3v.rds', compress = FALSE)

# Display script execution time and completion message
end_time <- Sys.time()
cat("Script completed at: ", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Total execution time: ", round(difftime(end_time, start_time, units = "mins"), 2), " minutes\n")
