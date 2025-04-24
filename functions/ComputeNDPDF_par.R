#' Compute N-dimensional PDFs for Two Time Periods
#'
#' This function computes n-dimensional Probability Density Functions (PDFs)
#' for a set of climate variables across multiple climate models. The PDFs are
#' computed for two time periods (present and future) and for each grid point
#' based on the provided longitude and latitude coordinates.
#'
#' @param variables A character vector of climate variables (e.g., 'pr', 'tas', 'psl').
#' @param model_names A character vector of climate model names to be processed.
#' @param data_dir A character string specifying the directory where the climate model data is stored.
#' @param year_present A numeric vector specifying the years representing the present time period.
#' @param year_future A numeric vector specifying the years representing the future time period.
#' @param lon A numeric vector of longitude coordinates.
#' @param lat A numeric vector of latitude coordinates.
#' @param range_var An array of dimension `[lon, lat, nvar, 2]` containing the range (min/max) for each variable.
#' @param nbins A numeric value specifying the number of bins to use when computing histograms for each variable.
#'
#' @return A list containing two arrays:
#' \item{present}{The n-dimensional PDF array for the present period.}
#' \item{future}{The n-dimensional PDF array for the future period.}
#' The arrays have dimensions `[lon, lat, nbins^nvar, num_models]`, where:
#' - `lon` and `lat` are the grid points.
#' - `nbins^nvar` is the number of bins for the joint PDF of all variables.
#' - `num_models` is the number of climate models.
#'
#' @examples
#' # Define parameters
#' variables <- c('pr', 'tas', 'psl')
#' model_names <- c('model1', 'model2', 'model3')
#' data_dir <- 'data/CMIP6_merged_all/'
#' year_present <- 1960:1990
#' year_future <- 2070:2100
#' lon <- 0:359
#' lat <- -90:90
#' range_var <- array(NA, dim = c(length(lon), length(lat), length(variables), 2))  # Example range array
#' nbins <- 20
#'
#' # Compute n-dimensional PDFs for both time periods
#' pdf_result <- compute_nd_pdf_optimized(variables, model_names, data_dir, year_present, year_future, lon, lat, range_var, nbins)
#'
#' @import ncdf4 future future.apply
#' @export
compute_nd_pdf_optimized <- function(variables, model_names, data_dir, year_present, year_future, lon, lat, range_var, nbins, workers) {
  n_var <- length(variables)  # Number of variables
  num_models <- length(model_names)
  pdf_matrix_present <- array(NA, dim = c(length(lon), length(lat), nbins^n_var, num_models))  # PDF matrix for present period
  pdf_matrix_future <- array(NA, dim = c(length(lon), length(lat), nbins^n_var, num_models))   # PDF matrix for future period

  # Set up parallel processing for each model
  plan(multisession, workers = workers)  # Use a limited number of workers
  options(future.globals.maxSize = 8 * 1024^3)  # Allow up to 8 GiB for exporting globals
  plan(sequential)  # Reset to sequential

  # Process models in parallel for each time period
  pdf_matrix_list <- future_lapply(seq_along(model_names), function(m) {
    model_name <- model_names[m]
    cat(paste0("Processing model: ", model_name, " at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))

    var_data_list_present <- list()
    var_data_list_future <- list()

    for (v in seq_along(variables)) {
      # Construct the file path
      file_path <- paste0(data_dir, model_name, '/', variables[v], '/', list.files(path = paste0(data_dir, model_name, '/', variables[v], '/'), pattern = glob2rx(paste0(variables[v], "_", model_name, "*.nc")))[1])
      nc_var <- nc_open(file_path)

      # Use the helper function to extract years from the time variable
      yyyy <- extract_years_from_time(nc_var)
      lon_var <- ncvar_get(nc_var, "lon")
      lat_var <- ncvar_get(nc_var, "lat")

      # Find indices for the longitude and latitude ranges
      lon_indices <- which(lon_var %in% lon)
      lat_indices <- which(lat_var %in% lat)
      start_lon <- min(lon_indices)
      start_lat <- min(lat_indices)

      # Get data slices for present and future periods
      iyyyy_present <- which(yyyy %in% year_present)
      var_data_present <- ncvar_get(nc_var, variables[v], start = c(start_lon, start_lat, min(iyyyy_present)), count = c(length(lon_indices), length(lat_indices), length(iyyyy_present)))

      iyyyy_future <- which(yyyy %in% year_future)
      var_data_future <- ncvar_get(nc_var, variables[v], start = c(start_lon, start_lat, min(iyyyy_future)), count = c(length(lon_indices), length(lat_indices), length(iyyyy_future)))

      nc_close(nc_var)

      # Apply log transform for precipitation if required
      if (variables[v] == 'pr') {
        var_data_present <- log(var_data_present + 1)
        var_data_future <- log(var_data_future + 1)
      }

      var_data_list_present[[v]] <- var_data_present
      var_data_list_future[[v]] <- var_data_future
    }

    # Compute PDFs for present period
    model_pdf_matrix_present <- array(NA, dim = c(length(lon), length(lat), nbins^n_var))
    for (i in seq_along(lon)) {
      for (j in seq_along(lat)) {
        pixel_data_present <- sapply(var_data_list_present, function(var) var[i, j, ])
        hist_tmp_present <- compute_histND(pixel_data_present, range_var[i, j, , ], nbins)
        model_pdf_matrix_present[i, j, ] <- hist_tmp_present / sum(hist_tmp_present)
      }
    }

    # Compute PDFs for future period
    model_pdf_matrix_future <- array(NA, dim = c(length(lon), length(lat), nbins^n_var))
    for (i in seq_along(lon)) {
      for (j in seq_along(lat)) {
        pixel_data_future <- sapply(var_data_list_future, function(var) var[i, j, ])
        hist_tmp_future <- compute_histND(pixel_data_future, range_var[i, j, , ], nbins)
        model_pdf_matrix_future[i, j, ] <- hist_tmp_future / sum(hist_tmp_future)
      }
    }

    return(list(present = model_pdf_matrix_present, future = model_pdf_matrix_future))
  })

  # Combine results into the main PDF matrices
  for (m in seq_along(model_names)) {
    pdf_matrix_present[, , , m] <- pdf_matrix_list[[m]]$present
    pdf_matrix_future[, , , m] <- pdf_matrix_list[[m]]$future
  }

  plan(sequential)  # Reset to sequential

  return(list(present = pdf_matrix_present, future = pdf_matrix_future))
}


# Helper function to extract years from time netcdf
extract_years_from_time <- function(nc_var) {
  # Determine the correct time dimension name
  time_dim <- if ("time" %in% names(nc_var$dim)) "time" else "valid_time"

  # Retrieve the raw time data
  time_raw <- ncvar_get(nc_var, time_dim)

  # Get the units of the time variable
  time_units <- ncatt_get(nc_var, time_dim, "units")$value

  # Check the time units format and convert accordingly
  if (grepl("since", time_units)) {
    reference_date <- as.Date(sub(".*since ", "", time_units))

    if (grepl("days", time_units)) {
      dates <- reference_date + time_raw
    } else if (grepl("seconds", time_units)) {
      dates <- reference_date + as.difftime(time_raw, units = "secs")
    }

    # Extract the year component
    return(as.numeric(format(dates, "%Y")))
  } else {
    stop("Unrecognized time units format.")
  }
}