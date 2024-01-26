OpenAndKDE1D <- function (model_names, variables,
                        year_present, year_future, period) {

  # Initialize data structures
  kde_matrix <- array(0, c(length(lon), length(lat), nbins1d,
                           length(model_names), length(variables)))

  # Loop through variables and models
  v <- 1
  for(var in variables){
    m <- 1
    for(model_name in model_names){
      dir_path <- paste0(data_dir, model_name, '/', var, '/')
      # Create the pattern
      pattern <- glob2rx(paste0(var, "_", model_name, "_", period, "*.nc"))

      # Get the filepath
      file_name <- list.files(path = dir_path,
                              pattern = pattern)
      file_path <- paste0(dir_path, file_name)
      print(file_path)

      # Check that there is only one matching file
      if (length(file_path) == 1) {
        nc <<- nc_open(file_path)
        # Extract and average data for present
        yyyy <- substr(as.character(nc.get.time.series(nc)), 1, 4)
        iyyyy <- which(yyyy %in% year_present)


        # For each grid point...
        for (i in seq_along(lon)){
          for (j in seq_along(lat)){
            print(c(i, j, var))
            tmp <- ncvar_get(nc, var, start = c(i, j, min(iyyyy)),
                             count = c(1, 1, length(iyyyy)))

            if(var == 'pr' ) {
              # Changing units to mm/day and log transform
              tmp <- log10((tmp * 86400) + 1)
            }

            # For the first model
            if(m == 1 ) {
              # Compute the range of the variable for grid-point i-j
              range_var <- range(tmp)
              # Calculate the 10% margin
              margin <- diff(range_var) * 0.1

              # Set the lower and upper bounds with a 10% margin
              range_var[1] <- range_var[1] - margin
              range_var[2] <- range_var[2] + margin
            }

            # TODO should we fix the bandwidth?
            dens_tmp <- density(tmp,
                                from = range_var[1], to = range_var[2],
                                n = nbins1d)

            kde_matrix[i, j, , m, v] <- dens_tmp$y

          }
        }

        # Create dimensions and initialize matrices for present and future data if this is the first model and variable
        if(j == 1 && i == 1) {
          lat <- ncvar_get(nc, "lat")
          lon <- ncvar_get(nc, "lon")
          data_matrix$present <- data_matrix$future <-
            array(0, c(length(lon), length(lat),
                       length(model_names), length(variables)))
        }

        # Extract and average data for present
        yyyy <- substr(as.character(nc.get.time.series(nc)), 1, 4)
        iyyyy <- which(yyyy %in% year_present)
        tmp <- apply(
          ncvar_get(nc, var, start = c(1, 1, min(iyyyy)), count = c(-1, -1, length(iyyyy))),
          1:2,
          mean
        )
        data_matrix$present[, , i, j] <- tmp
        data[['present']][[var]][[model_name]] <- tmp

        # Extract and average data for future
        iyyyy <- which(yyyy %in% year_future)
        tmp <- apply(
          ncvar_get(nc, var, start = c(1, 1, min(iyyyy)), count = c(-1, -1, length(iyyyy))),
          1:2,
          mean
        )
        data_matrix$future [, , i, j] <- tmp
        data[['future']][[var]][[model_name]] <- tmp

        # Close the file
        nc_close(nc)

      }else {
        # Handle the case where there are multiple or no matching files
        print("Error: Found multiple or no matching files")
      }
      # Update counter for models
      m <- m + 1
    }
    # Update counter for variables
    v <- v + 1
  }


  # Remove the counters
  remove(m, v)

  # Set the dimnames of the matrices
  dimnames(data_matrix$present) <- dimnames(data_matrix$future) <-
    list(lon = lon, lat = lat, model = paste0(model_names, model_names))

  # Return the output as a list of two elements
  output <- list(data, data_matrix, kde_matrix)
  return(output)
}