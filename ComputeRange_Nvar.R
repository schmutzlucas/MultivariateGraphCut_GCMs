# -----------------------------------------------
# Generic Script for Calculating and Merging Variable Ranges
# across Multiple Climate Models Using Parallel Processing
# -----------------------------------------------

# This script performs the following operations:
# 1. Installs and loads necessary R packages.
# 2. Loads custom functions from local directories.
# 3. Defines global variables such as longitude, latitude, and models.
# 4. Sets up parallel backend to speed up processing.
# 5. Defines a function to compute the range (min and max) values
#    for each variable across a list of climate models.
# 6. Uses parallel processing (`future_lapply`) to compute ranges.
# 7. Merges the computed ranges into a final global range for each variable.
# 8. Saves the results as RDS files.
# 9. Displays the total execution time.
# -----------------------------------------------

# --------- 1. Install and Load Necessary Libraries ---------
# Read the required package list from a file and install if not available
required_packages <- read.table("package_list.txt", sep = "\n")$V1
missing_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(missing_packages)) install.packages(missing_packages, repos = "https://cloud.r-project.org")

# Load required libraries
invisible(lapply(required_packages, library, character.only = TRUE))

# --------- 2. Load Custom Functions from Local Directory ---------
source_code_dir <- 'functions/'  # Define function directory path
invisible(lapply(list.files(source_code_dir, full.names = TRUE), source))

# --------- 3. Global Variables Definition ---------
# Define constants and data paths
lon <- 0:359                     # Longitude range
lat <- -70:70                    # Latitude range
year_interest <- 1960:2022       # Years of interest for analysis
data_dir <- 'data/CMIP6_merged_all/'  # Directory containing climate model data

# Automatically load variables and models from files
variables <- c('pr', 'tas', 'psl')   # List of variables dynamically loaded
model_names <- read.table('model_names_pr_tas_psl.txt')$V1  # List of models dynamically loaded

# --------- 4. Set Up Parallel Backend Using `future` ---------
plan(multisession, workers = 8)  # Set up parallel processing with 8 workers

# --------- 5. Function to Calculate Min and Max Ranges for Any Variable ---------
calculate_ranges <- function(variable, models, data_dir, years, lon_grid, lat_grid) {
  # Create an empty array to store min and max ranges for each model
  range_array <- array(NA, dim = c(length(lon_grid), length(lat_grid), 2, length(models)))

  for (m in seq_along(models)) {
    model_name <- models[m]
    # Construct the file path using the variable and model name dynamically
    file_pattern <- paste0(variable, "_", model_name, "*.nc")
    file_path <- paste0(data_dir, model_name, '/', variable, '/', list.files(paste0(data_dir, model_name, '/', variable, '/'), pattern = glob2rx(file_pattern))[1])

    # Open the NetCDF file and extract data
    nc_data <- nc_open(file_path)
    time_series <- substr(as.character(nc.get.time.series(nc_data)), 1, 4)
    year_indices <- which(time_series %in% years)

    # Extract spatial indices and read variable data
    lon_indices <- which(ncvar_get(nc_data, "lon") %in% lon_grid)
    lat_indices <- which(ncvar_get(nc_data, "lat") %in% lat_grid)

    # Retrieve variable data within the defined spatial and temporal bounds
    grid_data <- ncvar_get(nc_data, variable, start = c(min(lon_indices), min(lat_indices), min(year_indices)),
                           count = c(length(lon_indices), length(lat_indices), length(year_indices)))

    # Apply log transformation for precipitation variable
    if (variable == 'pr') grid_data <- log(grid_data + 1)

    # Calculate and store min and max values for each grid cell
    range_array[ , , 1, m] <- apply(grid_data, c(1, 2), min, na.rm = TRUE)
    range_array[ , , 2, m] <- apply(grid_data, c(1, 2), max, na.rm = TRUE)

    nc_close(nc_data)  # Close NetCDF file after processing
  }
  return(range_array)
}

# --------- 6. Calculate Ranges for Each Variable in Parallel ---------
range_results <- future_lapply(variables, calculate_ranges, models = model_names, data_dir = data_dir,
                               years = year_interest, lon_grid = lon, lat_grid = lat)

# --------- 7. Assign Results Dynamically to Named List ---------
# Create a named list to store range arrays for each variable
range_vars <- setNames(range_results, variables)

# --------- 8. Save Each Variable's Range Separately ---------
# Save ranges for each variable in RDS files using dynamic names
lapply(variables, function(var) saveRDS(range_vars[[var]], paste0('ranges/', var, '_range_AllModelsPar_1950-2022_70deg.rds')))

# --------- 9. Merge Results into Final Global Ranges ---------
# Initialize final range array for each variable
final_ranges <- lapply(variables, function(x) array(NA, dim = c(length(lon), length(lat), 2)))

# Compute global min and max for each grid cell across all models
for (i in seq_along(lon)) {
  for (j in seq_along(lat)) {
    for (v in seq_along(variables)) {
      final_ranges[[v]][i, j, ] <- range(range_vars[[variables[v]]][i, j, , ], na.rm = TRUE)
    }
  }
}

# --------- 10. Save the Final Merged Ranges ---------
# Save final ranges as a single RDS file
saveRDS(final_ranges, 'ranges/range_var_final_allModelsPar_1950-2022_70deg.rds')

# --------- 11. Display Script Execution Time ---------
end_time <- Sys.time()
cat("Script completed at: ", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Total execution time: ", round(difftime(end_time, start_time, units = "mins"), 2), " minutes\n")
