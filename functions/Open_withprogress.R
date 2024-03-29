#' Open and compute 1D KDE of climate model data
#'
#' This function opens NetCDF climate model data files and computes the 1D Kernel Density Estimation (KDE) for each grid point and specified variables.
#'
#' @param model_names A character vector containing the names of the climate models to be processed.
#' @param variables A character vector containing the names of the variables to be processed.
#' @param year_present A character vector containing the years to be considered for the present period.
#' @param year_future A character vector containing the years to be considered for the future period.
#' @param period A character string specifying the period of the data (e.g., "historical" or "rcp85").
#'
#' @return A list containing the computed 1D KDE matrix and the variable ranges.
#'
#' @examples
#' # Example usage:
#' # model_names <- c("Model1", "Model2")
#' # variables <- c("pr", "tas")
#' # year_present <- 1979:1998
#' # year_future <- 1999:2019
#' # period <- "historical"
#' # results <- OpenAndKDE1D_new(model_names, variables, year_present, year_future, period)
#'
#' @import ncdf4
#' @export
OpenAndKDE1D_new <- function (model_names, variables,
                              year_present, year_future, period) {

  # Initialize data structures
  pdf_matrix <- array(0, c(length(lon), length(lat), nbins1d,
                           length(model_names), length(variables)))
  # Initialize data structures
  mids_matrix <- array(0, c(length(lon), length(lat), nbins1d,
                            length(model_names), length(variables)))

  range_var <- list()



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

        # Get the entire 2D-time model as array
        tmp_grid <- ncvar_get(nc, var, start = c(1, 1, min(iyyyy)), count = c(-1, -1, length(iyyyy)))

        # For each grid point...
        for (i in seq_along(lon)) {
          for (j in seq_along(lat)) {
            if (i%%10 == 0 && j %% 100 == 0) {
              print(c(v, var, m, model_name, i, j))
            }
            tmp <- tmp_grid[i, j, ]

            if (var == 'pr') {
              tmp <- log2((tmp * 86400) + 1)
            }

            if (m == 1) {
              if (i == 1 && j == 1) {
                range_var[[var]] <- array(data = NA, dim = c(length(lon), length(lat), 2))
              }
              if (var == 'pr') {
                range_var[[var]][i, j, 1] <- 0
                range_var[[var]][i, j, 2] <- range(tmp)[2] * 2
              }
              else{
                range_var[[var]][i, j, 1] <- range(tmp)[1] - diff(range(tmp)) * 0.3
                range_var[[var]][i, j, 2] <- range(tmp)[2] + diff(range(tmp)) * 0.3
              }
            }
            # Compute the breaks
            breaks <- seq(from = range_var[[var]][i, j, 1],
                          to = range_var[[var]][i, j, 2],
                          length.out = nbins1d + 1)

            # Replace values lower or higher than the range with min or max
            min_range <- range_var[[var]][i, j, 1]
            max_range <- range_var[[var]][i, j, 2]
            tmp[tmp < min_range] <- min_range
            tmp[tmp > max_range] <- max_range

            # Compute the histogram using the modified data
            dens_tmp <- hist(tmp, breaks = breaks, plot = FALSE)

            pdf_matrix[i, j, , m, v] <- dens_tmp$counts / sum(dens_tmp$counts)
            mids_matrix[i,j,,m,v] <- dens_tmp$mids

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