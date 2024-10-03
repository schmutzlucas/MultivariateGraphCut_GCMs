# Install and load necessary packages
list_of_packages <- readLines("package_list.txt")

# Install missing packages
new_packages <- setdiff(list_of_packages, rownames(installed.packages()))
if (length(new_packages)) {
  install.packages(new_packages, repos = "https://cloud.r-project.org")
}

# Load all packages
lapply(list_of_packages, library, character.only = TRUE)

# Conditionally install and load 'gcoWrapR'
if (!requireNamespace("gcoWrapR", quietly = TRUE)) {
  devtools::install_github("thaos/gcoWrapR")
}
library(gcoWrapR)

# Load additional packages used in the script
library(future)
library(future.apply)
library(progressr)
library(logger)

# Define the function to calculate ranges
calculate_ranges <- function(variable, model_names, data_dir, year_interest, lon, lat) {
  num_models <- length(model_names)
  lon_len <- length(lon)
  lat_len <- length(lat)
  range_var <- array(NA_real_, dim = c(lon_len, lat_len, 2, num_models))

  for (m in seq_len(num_models)) {
    model_name <- model_names[m]
    file_dir <- file.path(data_dir, model_name, variable)
    file_pattern <- paste0(variable, "_", model_name, "*.nc")
    file_list <- list.files(path = file_dir, pattern = glob2rx(file_pattern), full.names = TRUE)

    if (length(file_list) == 0) {
      log_warn("No files found for variable '{variable}' and model '{model_name}'. Skipping.")
      next
    }

    file_path <- file_list[1]  # Modify as needed

    nc_var <- tryCatch(nc_open(file_path), error = function(e) {
      log_error("Error opening NetCDF file '{file_path}': {e$message}")
      return(NULL)
    })
    if (is.null(nc_var)) next

    # Extract time variable
    time_var <- ncvar_get(nc_var, "time")
    time_units <- ncatt_get(nc_var, "time", "units")$value
    time_calendar <- ncatt_get(nc_var, "time", "calendar")$value

    time_dates <- nc.get.time.series(nc_var, v = "time", time.dim.name = "time")
    yyyy <- format(time_dates, "%Y")
    iyyyy <- which(yyyy %in% year_interest)
    if (length(iyyyy) == 0) {
      log_warn("No matching years found in file '{file_path}'. Skipping.")
      nc_close(nc_var)
      next
    }

    # Extract spatial variables
    lon_var <- ncvar_get(nc_var, "lon")
    lat_var <- ncvar_get(nc_var, "lat")

    lon_indices <- which(lon_var %in% lon)
    lat_indices <- which(lat_var %in% lat)
    if (length(lon_indices) == 0 || length(lat_indices) == 0) {
      log_warn("No matching lon/lat indices in file '{file_path}'. Skipping.")
      nc_close(nc_var)
      next
    }

    # Read the data
    start <- c(min(lon_indices), min(lat_indices), min(iyyyy))
    count <- c(length(lon_indices), length(lat_indices), length(iyyyy))
    tmp_grid_var <- ncvar_get(nc_var, variable, start = start, count = count)
    nc_close(nc_var)

    if (variable == 'pr') tmp_grid_var <- log(tmp_grid_var + 1)

    # Compute min and max across the time dimension
    range_var[ , , 1, m] <- apply(tmp_grid_var, c(1, 2), min, na.rm = TRUE)
    range_var[ , , 2, m] <- apply(tmp_grid_var, c(1, 2), max, na.rm = TRUE)
  }

  return(range_var)
}

# Main script execution
main <- function() {
  # Set up logging
  log_appender(appender_console())
  start_time <- Sys.time()
  log_info("Script started at: {format(start_time, '%Y-%m-%d %H:%M:%S')}")

  # Loading local functions
  source_code_dir <- 'functions/'
  file_paths <- list.files(source_code_dir, pattern = "\\.R$", full.names = TRUE)
  lapply(file_paths, source)

  # Global variables
  lon <- 0:359
  lat <- -70:70
  year_interest <- 1950:2022
  data_dir <- 'data/CMIP6_merged_all/'
  variables <- c('pr', 'tas', 'psl')
  model_names <- read.table('model_names_pr_tas_psl.txt', stringsAsFactors = FALSE)$V1

  # Set up parallel backend using 'future'
  plan(multisession)

  # Set up progress bar
  handlers(global = TRUE)

  # Calculate ranges for each variable asynchronously
  with_progress({
    p <- progressor(along = variables)
    range_results <- future_lapply(variables, function(var) {
      p(sprintf("Processing variable: %s", var))
      calculate_ranges(var, model_names, data_dir, year_interest, lon, lat)
    })
  })

  # Assign names to range_results
  names(range_results) <- variables

  # Save each variable's range separately
  for (var in variables) {
    dir.create("ranges", showWarnings = FALSE)
    file_name <- sprintf('ranges/%s_range_AllModelsPar_%d-%d_70deg.rds', var, min(year_interest), max(year_interest))
    saveRDS(range_results[[var]], file_name)
  }

  # Final global range computation by merging results
  range_var_final <- list()
  for (v in seq_along(variables)) {
    var <- variables[v]
    range_min <- apply(range_results[[var]][ , , 1, , drop = FALSE], c(1, 2), min, na.rm = TRUE)
    range_max <- apply(range_results[[var]][ , , 2, , drop = FALSE], c(1, 2), max, na.rm = TRUE)
    range_var_final[[var]] <- array(NA_real_, dim = c(length(lon), length(lat), 2))
    range_var_final[[var]][ , , 1] <- range_min
    range_var_final[[var]][ , , 2] <- range_max
  }

  # Save the final merged ranges
  final_file_name <- sprintf('ranges/range_var_final_allModelsPar_%d-%d_70deg.rds', min(year_interest), max(year_interest))
  saveRDS(range_var_final, final_file_name)

  # Log script completion
  end_time <- Sys.time()
  log_info("Script completed at: {format(end_time, '%Y-%m-%d %H:%M:%S')}")
  log_info("Total execution time: {round(difftime(end_time, start_time, units = 'mins'), 2)} minutes")
}

# Call the main function
main()
