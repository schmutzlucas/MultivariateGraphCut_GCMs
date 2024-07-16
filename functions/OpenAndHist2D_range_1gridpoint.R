OpenAndHist2D_range_lonlat <- function (model_names, variables, year_interest, range_var, lon_vec, lat_vec) {
  # Initialize data structures
  pdf_matrix <- array(0, c(length(lon_vec), length(lat_vec), nbins1d^2, length(model_names)))
  x_breaks <- array(0, c(length(lon_vec), length(lat_vec), nbins1d+1, length(model_names)))
  y_breaks <- array(0, c(length(lon_vec), length(lat_vec), nbins1d+1, length(model_names)))

  # Loop through variables and models
  m <- 1
  for(model_name in model_names){
    # Variable 1
    dir_path <- paste0(data_dir, model_name, '/', variables[1], '/')
    # Create the pattern
    pattern <- glob2rx(paste0(variables[1], "_", model_name, "*.nc"))

    # Get the filepath
    file_name <- list.files(path = dir_path, pattern = pattern, full.names = TRUE)
    stopifnot(length(file_name) == 1)
    nc1 <- nc_open(file_name)

    # Variable 2
    dir_path <- paste0(data_dir, model_name, '/', variables[2], '/')
    # Create the pattern
    pattern <- glob2rx(paste0(variables[2], "_", model_name, "*.nc"))

    # Get the filepath
    file_name <- list.files(path = dir_path, pattern = pattern, full.names = TRUE)
    stopifnot(length(file_name) == 1)
    nc2 <- nc_open(file_name)

    # Extract and average data for present
    yyyy2 <- substr(as.character(nc.get.time.series(nc2)), 1, 4)
    iyyyy2 <- which(yyyy2 %in% year_interest)

    # Get the data only for the specified lon and lat points
    i <- 1
    j <- 1
    print(file_name)
    print(length(lon_vec))
    for (k in 1:length(lon_vec)) {
      print(lon_vec[k])
      cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")


        tmp1 <- ncvar_get(nc1, variables[1], start = c(lon_vec[k], lat_vec[k], min(iyyyy2)), count = c(1, 1, length(iyyyy2)))
        tmp2 <- ncvar_get(nc2, variables[2], start = c(lon_vec[k], lat_vec[k], min(iyyyy2)), count = c(1, 1, length(iyyyy2)))

        if (variables[1] == 'pr') {
          tmp1 <- log(tmp1 + 1)
        }
        if (variables[2] == 'pr') {
          tmp2 <- log(tmp2 + 1)
        }

        # Compute the breaks
        breaks1 <- seq(from = range_var[[variables[1]]][lon_vec[k], lat_vec[k], 1], to = range_var[[variables[1]]][lon_vec[k], lat_vec[k], 2], length.out = nbins1d + 1)
        breaks2 <- seq(from = range_var[[variables[2]]][lon_vec[k], lat_vec[k], 1], to = range_var[[variables[2]]][lon_vec[k], lat_vec[k], 2], length.out = nbins1d + 1)

        # Compute the histogram using the modified data
        dens_tmp <- hist2(tmp1, tmp2, xbreaks = breaks1, ybreaks = breaks2, plot = FALSE)
        dens_tmp$z <- replace(dens_tmp$z, is.nan(dens_tmp$z), 0)
        hist_tmp <- dens_tmp$z
        hist_tmp <- replace(hist_tmp, is.na(hist_tmp), 0)
        pdf_matrix[i, j, , m] <- c(hist_tmp / sum(hist_tmp))
        x_breaks[i,j,,m] <- dens_tmp$x
        y_breaks[i,j,,m] <- dens_tmp$y
      i <- i + 1
      j <- j + 1
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