OpenAndHist2D_range_lonlat_2 <- function (model_names, variables, year_interest, range_var, lon_vec, lat_vec) {
  # Initialize data structures
  pdf_matrix <- array(0, c(length(lon_vec), length(lat_vec), nbins1d^2, length(model_names)))
  x_breaks <- array(0, c(length(lon_vec), length(lat_vec), nbins1d+1, length(model_names)))
  y_breaks <- array(0, c(length(lon_vec), length(lat_vec), nbins1d+1, length(model_names)))

  # Loop through variables and models
  m <- 1
  for(model_name in model_names) {
    # Variable 1
    dir_path <- paste0(data_dir, model_name, '/', variables[1], '/')
    pattern <- glob2rx(paste0(variables[1], "_", model_name, "*.nc"))
    file_name <- list.files(path = dir_path, pattern = pattern, full.names = TRUE)
    stopifnot(length(file_name) == 1)
    nc1 <- nc_open(file_name)

    # Variable 2
    dir_path <- paste0(data_dir, model_name, '/', variables[2], '/')
    pattern <- glob2rx(paste0(variables[2], "_", model_name, "*.nc"))
    file_name <- list.files(path = dir_path, pattern = pattern, full.names = TRUE)
    stopifnot(length(file_name) == 1)
    nc2 <- nc_open(file_name)

    # Extract and average data for present
    yyyy2 <- substr(as.character(nc.get.time.series(nc2)), 1, 4)
    iyyyy2 <- which(yyyy2 %in% year_interest)

    # Get variable dimensions
    lon_dim1 <- length(ncvar_get(nc1, "lon"))
    lat_dim1 <- length(ncvar_get(nc1, "lat"))
    lon_dim2 <- length(ncvar_get(nc2, "lon"))
    lat_dim2 <- length(ncvar_get(nc2, "lat"))


    print(file_name)
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    for (i in 1:length(lon_vec)) {
      for (j in 1:length(lat_vec)) {

        # Validate indices
        if (lon_vec[i] > lon_dim1 || lat_vec[j] > lat_dim1 || lon_vec[i] < 1 || lat_vec[j] < 1) {
          stop("Error: Index exceeds dimension bound for lon/lat")
        }

        # Get the data only for the specified lon and lat points
        tmp1 <- ncvar_get(nc1, variables[1], start = c(lon_vec[i], lat_vec[j], min(iyyyy2)), count = c(1, 1, length(iyyyy2)))
        tmp2 <- ncvar_get(nc2, variables[2], start = c(lon_vec[i], lat_vec[j], min(iyyyy2)), count = c(1, 1, length(iyyyy2)))

        if (variables[1] == 'pr') {
          tmp1 <- log(tmp1 + 1)
        }
        if (variables[2] == 'pr') {
          tmp2 <- log(tmp2 + 1)
        }

        # Compute the breaks
        breaks1 <- seq(from = range_var[[variables[1]]][lon_vec[i], lat_vec[j], 1], to = range_var[[variables[1]]][lon_vec[i], lat_vec[j], 2], length.out = nbins1d + 1)
        breaks2 <- seq(from = range_var[[variables[2]]][lon_vec[i], lat_vec[j], 1], to = range_var[[variables[2]]][lon_vec[i], lat_vec[j], 2], length.out = nbins1d + 1)

        # Compute the histogram using the modified data
        dens_tmp <- hist2(tmp1, tmp2, xbreaks = breaks1, ybreaks = breaks2, plot = FALSE)
        dens_tmp$z <- replace(dens_tmp$z, is.nan(dens_tmp$z), 0)
        hist_tmp <- dens_tmp$z
        hist_tmp <- replace(hist_tmp, is.na(hist_tmp), 0)
        pdf_matrix[i, j, , m] <- c(hist_tmp / sum(hist_tmp))
        x_breaks[i,j,,m] <- dens_tmp$x
        y_breaks[i,j,,m] <- dens_tmp$y
      }
    }

    # Close the file
    nc_close(nc1)
    nc_close(nc2)
    # Update counter for models
    m <- m + 1
  }

  # Remove the counters
  remove(m)

  # Return the output as a list of two elements
  output <- list(pdf_matrix, range_var, x_breaks, y_breaks, tmp1, tmp2)
  return(output)
}
