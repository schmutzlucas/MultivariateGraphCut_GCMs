OpenAndKDE1D_new <- function (model_names, variables,
                              year_present, year_future, period) {

  # Initialize data structures
  # kde_matrix <- array(0, c(length(lon), length(lat), nbins1d,
  #                          length(model_names), length(variables)))
  # Initialize data structures
  # range_var <- list()


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

        # Get the entire 3D matrix
        tmp_grid <- ncvar_get(nc, var, start = c(1, 1, min(iyyyy)), count = c(-1, -1, length(iyyyy)))

        # For each grid point...
        for (i in seq_along(lon)) {
          for (j in seq_along(lat)) {
            if (i%%10 == 0 && j %% 100 == 0) {
              print(c(v, var, m, model_name, i, j))
            }
            tmp <- tmp_grid[i, j, ]

            if (var == 'pr') {
              tmp <- log10((tmp * 86400) + 1)
            }

            if (m == 1) {
              if (i == 1 && j == 1) {
                range_var[[var]] <<- array(data = NA, dim = c(length(lon), length(lat), 2))
              }
              range_var[[var]][i, j, 1] <<- range(tmp)[1] - diff(range(tmp)) * 0.1
              range_var[[var]][i, j, 2] <<- range(tmp)[2] + diff(range(tmp)) * 0.1
            }

            dens_tmp <- density(tmp,
                                from = range_var[[var]][i, j, 1],
                                to = range_var[[var]][i, j, 2],
                                n = nbins1d)

            pdf_matrix[i, j, , m, v] <<- dens_tmp$y * dens_tmp$bw
          }
        }

        # # Create dimensions and initialize matrices for present and future data if this is the first model and variable
        # if(j == 1 && i == 1) {
        #   lat <- ncvar_get(nc, "lat")
        #   lon <- ncvar_get(nc, "lon")
        #   data_matrix$present <- data_matrix$future <-
        #     array(0, c(length(lon), length(lat),
        #                length(model_names), length(variables)))
        # }
        #
        # # Extract and average data for present
        # yyyy <- substr(as.character(nc.get.time.series(nc)), 1, 4)
        # iyyyy <- which(yyyy %in% year_present)
        # tmp <- apply(
        #   ncvar_get(nc, var, start = c(1, 1, min(iyyyy)), count = c(-1, -1, length(iyyyy))),
        #   1:2,
        #   mean
        # )
        # data_matrix$present[, , i, j] <- tmp
        # data[['present']][[var]][[model_name]] <- tmp
        #
        # # Extract and average data for future
        # iyyyy <- which(yyyy %in% year_future)
        # tmp <- apply(
        #   ncvar_get(nc, var, start = c(1, 1, min(iyyyy)), count = c(-1, -1, length(iyyyy))),
        #   1:2,
        #   mean
        # )
        # data_matrix$future [, , i, j] <- tmp
        # data[['future']][[var]][[model_name]] <- tmp

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

  # Return the output as a list of two elements
  output <- list(pdf_matrix, range_var)
  return(output)
}