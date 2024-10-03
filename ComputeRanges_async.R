# Install and load necessary libraries
list_of_packages <- read.table("package_list.txt", sep = "\n")$V1
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages, repos = "https://cloud.r-project.org")

library(future)
library(future.apply)
library(progressr)
library(ncdf4)
library(devtools)
lapply(list_of_packages, library, character.only = TRUE)
install_github("thaos/gcoWrapR")

# Set up `progressr` for better visibility of parallel progress
handlers(global = TRUE)  # Enable global progress handlers
handlers("txtprogressbar")  # Use a text progress bar for console updates

# Loading local functions
source_code_dir <- 'functions/'  # The directory where all functions are saved.
file_paths <- list.files(source_code_dir, full.names = TRUE)
for (path in file_paths) {
  source(path)
}

# Global variables
lon <- 0:359
lat <- -70:70
year_interest <- 1950:2022
data_dir <- 'data/CMIP6_merged_all/'
variables <- c('pr', 'tas', 'psl')
model_names <- read.table('model_names_pr_tas_psl.txt')$V1

# Set up parallel backend using `future`
plan(multisession, workers = 3)  # Use 3 separate R sessions

# Function for calculating ranges for each variable across all models
calculate_ranges <- function(variable, model_names, data_dir, year_interest, lon, lat) {
  # Start progress tracking
  p <- progressr::progressor(steps = length(model_names))
  message(paste0("Starting range calculation for variable: '", variable, "' at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

  range_var <- array(data = NA, dim = c(length(lon), length(lat), 2, length(model_names)))  # Initialize result array

  for (m in seq_along(model_names)) {
    model_name <- model_names[m]
    dir_path <- paste0(data_dir, model_name, '/', variable, '/')
    pattern <- glob2rx(paste0(variable, "_", model_name, "*.nc"))
    file_name <- list.files(path = dir_path, pattern = pattern)
    file_path <- paste0(dir_path, file_name)

    # Log model and file information with progress
    p(sprintf("Processing Model '%s' for Variable '%s' [File: %s]", model_name, variable, basename(file_path)))
    message(paste0("[INFO] Processing Model: '", model_name, "' for Variable: '", variable, "' at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    message(paste0("[INFO] File Path: ", file_path))

    # Open the NetCDF file
    nc_var <- nc_open(file_path)

    # Extract time indices and data
    yyyy <- substr(as.character(nc.get.time.series(nc_var)), 1, 4)
    iyyyy <- which(yyyy %in% year_interest)
    lon_var <- ncvar_get(nc_var, "lon")
    lat_var <- ncvar_get(nc_var, "lat")

    # Match longitude and latitude indices
    lon_indices <- which(lon_var %in% lon)
    lat_indices <- which(lat_var %in% lat)
    start_lon <- min(lon_indices)
    count_lon <- length(lon_indices)
    start_lat <- min(lat_indices)
    count_lat <- length(lat_indices)

    # Extract data slice for the specified variable
    tmp_grid_var <- ncvar_get(nc_var, variable, start = c(start_lon, start_lat, min(iyyyy)),
                              count = c(count_lon, count_lat, length(iyyyy)))

    # Apply log transformation for precipitation variable
    if (variable == 'pr') {
      tmp_grid_var <- log(tmp_grid_var + 1)
    }

    # Calculate min and max for each pixel
    range_var[ , , 1, m] <- apply(tmp_grid_var, c(1, 2), min, na.rm = TRUE)
    range_var[ , , 2, m] <- apply(tmp_grid_var, c(1, 2), max, na.rm = TRUE)

    # Close the NetCDF file
    nc_close(nc_var)
    message(paste0("[INFO] Completed processing for Model: '", model_name, "' for Variable: '", variable, "' at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    p()
  }

  message(paste0("Completed range calculation for variable: '", variable, "' at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  return(range_var)
}

# Enhanced progress messages with `progressr`
message(paste0("Starting parallel range calculations for all variables at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

# Run `calculate_ranges` for each variable asynchronously using `progressr` wrapper
range_results <- with_progress({
  future_lapply(variables, function(var) {
    message(paste0("[TASK] Starting calculation for Variable: '", var, "' at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    result <- calculate_ranges(var, model_names, data_dir, year_interest, lon, lat)
    message(paste0("[TASK] Completed calculation for Variable: '", var, "' at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    return(result)
  })
})

# Remaining part of the code...


# Assign the results to respective variables
range_var_1 <- range_results[[1]]
range_var_2 <- range_results[[2]]
range_var_3 <- range_results[[3]]

# Save each variable's range separately
saveRDS(range_var_1, 'ranges/pr_log_range_AllModelsPar_1950-2022_70deg.rds', compress = FALSE)
saveRDS(range_var_2, 'ranges/tas_range_AllModelsPar_1950-2022_70deg.rds', compress = FALSE)
saveRDS(range_var_3, 'ranges/psl_range_AllModelsPar_1950-2022_70deg.rds', compress = FALSE)

message(paste0("Finished saving individual range results at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
flush.console()

# Final global range computation by merging results
range_var_final <- list()
for (var in variables) {
  range_var_final[[var]] <- array(data = NA, dim = c(length(lon), length(lat), 2))
}

# Compute final range by taking the min and max over all models
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

# Save the final merged ranges
saveRDS(range_var_final, 'ranges/range_var_final_allModelsPar_1950-2022_70deg.rds', compress = FALSE)

message(paste0("Completed full processing at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
flush.console()
