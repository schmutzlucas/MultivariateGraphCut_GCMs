# -----------------------------------------------
# Script to Calculate and Merge Variable Ranges
# across Multiple Climate Models (N-variable capable)
# -----------------------------------------------

# --------- 1. Install and Load Necessary Libraries ---------
required_packages <- read.table("package_list.txt", sep = "\n")$V1
missing_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(missing_packages)) install.packages(missing_packages, repos = "https://cloud.r-project.org")

# Load required libraries
invisible(lapply(required_packages, library, character.only = TRUE))

# Additional Libraries
library(future)
library(future.apply)
library(ncdf4)
library(devtools)
install_github("schmutzlucas/gcoWrapR")

# --------- 2. Load Custom Functions ---------
source_code_dir <- 'functions/'  # Define function directory path
invisible(lapply(list.files(source_code_dir, full.names = TRUE), source))

# Start timing the script execution
start_time <- Sys.time()
cat("Script started at: ", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")

# --------- 3. Global Variables ---------
lon <- 0:359  # Longitude range
lat <- -90:90  # Latitude range
year_interest <- 1950:2100  # Years of interest
data_dir <- 'data/CMIP6_summer_Apr15-Oct14/'  # Directory for climate data
variables <- c('pr', 'tas', 'psl')  # List of variables
model_names <- read.table('model_names_pr_tas_psl_perfect_model_without_duplicate.txt')$V1  # List of models

# --------- 4. Parallel Processing Setup ---------

n_cores <- 2

if (.Platform$OS.type == "unix") {
  # Linux and macOS
  plan(multicore, workers = n_cores)
} else {
  # Windows fallback
  plan(multisession, workers = n_cores)
}


# --------- 5. Function to Extract Years from Time NetCDF ---------
extract_years_from_time <- function(nc_var) {
  time_dim <- if ("time" %in% names(nc_var$dim)) "time" else "valid_time"
  time_raw <- ncvar_get(nc_var, time_dim)
  time_units <- ncatt_get(nc_var, time_dim, "units")$value

  if (grepl("since", time_units)) {
    reference_date <- as.Date(sub(".*since ", "", time_units))
    if (grepl("days", time_units)) {
      dates <- reference_date + time_raw
    } else if (grepl("seconds", time_units)) {
      dates <- reference_date + as.difftime(time_raw, units = "secs")
    }
    return(as.numeric(format(dates, "%Y")))
  } else {
    stop("Unrecognized time units format.")
  }
}

# --------- 6. Function to Compute Ranges for Each Variable ---------
calculate_ranges <- function(variable, models, data_dir, years, lon_grid, lat_grid) {
  range_array <- array(NA, dim = c(length(lon_grid), length(lat_grid), 2, length(models)))

  for (m in seq_along(models)) {
    model_name <- models[m]
    file_pattern <- paste0(variable, "_", model_name, "*.nc")
    file_path <- paste0(data_dir, model_name, '/', variable, '/', list.files(paste0(data_dir, model_name, '/', variable, '/'), pattern = glob2rx(file_pattern))[1])

    # Open NetCDF file and extract metadata
    nc_data <- nc_open(file_path)
    time_years <- extract_years_from_time(nc_data)

    # Filter the required year range
    year_indices <- which(time_years %in% years)

    # Extract longitude and latitude indices
    lon_indices <- which(ncvar_get(nc_data, "lon") %in% lon_grid)
    lat_indices <- which(ncvar_get(nc_data, "lat") %in% lat_grid)

    # Retrieve data slice
    tmp_grid_var <- ncvar_get(
      nc_data, variable,
      start = c(min(lon_indices), min(lat_indices), min(year_indices)),
      count = c(length(lon_indices), length(lat_indices), length(year_indices)),
      collapse_degen = FALSE
    )

    # Transform the variable data if necessary (e.g., log transformation for precipitation)
    if (variable == 'pr') tmp_grid_var <- log(tmp_grid_var + 1)

    # Use hyperslabs to directly calculate min and max for time dimension
    min_values <- apply(tmp_grid_var, c(1, 2), min, na.rm = TRUE)
    max_values <- apply(tmp_grid_var, c(1, 2), max, na.rm = TRUE)

    # Store the results in range_array
    range_array[ , , 1, m] <- min_values
    range_array[ , , 2, m] <- max_values

    nc_close(nc_data)  # Close NetCDF file
    gc()
  }

  return(range_array)
}

# --------- 7. Compute Ranges for All Variables ---------
range_var <- future_lapply(variables, calculate_ranges, models = model_names, data_dir = data_dir,
                           years = year_interest, lon_grid = lon, lat_grid = lat)

plan(sequential)

# --------- 8. Compute Final Global Ranges ---------
final_ranges <- lapply(seq_along(variables), function(v) {
  range_data <- range_var[[v]]
  global_range <- array(NA, dim = c(length(lon), length(lat), 2))

  for (i in seq_along(lon)) {
    for (j in seq_along(lat)) {
      global_range[i, j, 1] <- min(range_data[i, j, 1, ], na.rm = TRUE)  # Global min
      global_range[i, j, 2] <- max(range_data[i, j, 2, ], na.rm = TRUE)  # Global max
    }
  }
  return(global_range)
})

names(final_ranges) <- variables

# --------- 9. Save the Final Merged Ranges ---------
saveRDS(list( ranges = final_ranges), 'ranges/range_var_summer_CMIP6_1950-2100_3v.rds')

# --------- 10. Display Script Execution Time ---------
end_time <- Sys.time()
cat("Script completed at: ", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Total execution time: ", round(difftime(end_time, start_time, units = "mins"), 2), " minutes\n")
