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
OpenAndHist2D <- function (model_names, variables,
                              year_interest, period) {

  # Initialize data structures
  pdf_matrix <- array(0, c(length(lon), length(lat), nbins1d^2,
                           length(model_names)))
  x_breaks <- array(0, c(length(lon), length(lat), nbins1d+1,
                           length(model_names)))
  y_breaks <- array(0, c(length(lon), length(lat), nbins1d+1,
                           length(model_names)))
  range_var <- list()



  # Loop through variables and models
  m <- 1
  for(model_name in model_names){
    # Variable 1
    dir_path <- paste0(data_dir, model_name, '/', variables[1], '/')
    # Create the pattern
    pattern <- glob2rx(paste0(variables[1], "_", model_name, "_", period, "*.nc"))

    # Get the filepath
    file_name <- list.files(path = dir_path,
                            pattern = pattern)
    file_path <- paste0(dir_path, file_name)
    print(file_path)

    # Check that there is only one matching file
    nc1 <<- nc_open(file_path)

    # Extract and average data for present
    yyyy1 <- substr(as.character(nc.get.time.series(nc1)), 1, 4)
    iyyyy1 <- which(yyyy1 %in% year_interest)

    # Get the entire 2D-time model as array
    tmp_grid_1 <- ncvar_get(nc1, variables[1], start = c(1, 1, min(iyyyy1)), count = c(-1, -1, length(iyyyy1)))


    # Variable 2
    dir_path <- paste0(data_dir, model_name, '/', variables[2], '/')
    # Create the pattern
    pattern <- glob2rx(paste0(variables[2], "_", model_name, "_", period, "*.nc"))

    # Get the filepath
    file_name <- list.files(path = dir_path,
                            pattern = pattern)
    file_path <- paste0(dir_path, file_name)
    print(file_path)

    # Check that there is only one matching file
    nc2 <<- nc_open(file_path)

    # Extract and average data for present
    yyyy2 <- substr(as.character(nc.get.time.series(nc2)), 1, 4)
    iyyyy2 <- which(yyyy2 %in% year_interest)

    # Get the entire 2D-time model as array
    tmp_grid_2 <- ncvar_get(nc2, variables[2], start = c(1, 1, min(iyyyy2)), count = c(-1, -1, length(iyyyy2)))

    # For each grid point...
    for (i in seq_along(lon)) {
      for (j in seq_along(lat)) {
        if (i%%10 == 0 && j %% 100 == 0) {
          print(c(m, model_name, i, j))
        }
        tmp1 <- tmp_grid_1[i, j, ]
        tmp2 <- tmp_grid_2[i, j, ]

        if (variables[1] == 'pr') {
          tmp1 <- log2((tmp1 * 86400) + 1)
        }
        if (variables[2] == 'pr') {
          tmp2 <- log2((tmp2 * 86400) + 1)
        }

        if (m == 1) {

          if (i == 1 && j == 1) {
            range_var[[variables[1]]] <- array(data = NA, dim = c(length(lon), length(lat), 2))
            range_var[[variables[2]]] <- array(data = NA, dim = c(length(lon), length(lat), 2))
          }
          if (variables[1] == 'pr') {
            range_var[[variables[1]]][i, j, 1] <- 0
            range_var[[variables[1]]][i, j, 2] <- range(tmp1)[2] * 1.2

            range_var[[variables[2]]][i, j, 1] <- range(tmp2)[1] - diff(range(tmp2)) * 0.3
            range_var[[variables[2]]][i, j, 2] <- range(tmp2)[2] + diff(range(tmp2)) * 0.3
          }
          else if (variables[2] == 'pr') {
            range_var[[variables[2]]][i, j, 1] <- 0
            range_var[[variables[2]]][i, j, 2] <- range(tmp2)[2] * 1.2

            range_var[[variables[1]]][i, j, 1] <- range(tmp1)[1] - diff(range(tmp1)) * 0.3
            range_var[[variables[1]]][i, j, 2] <- range(tmp1)[2] + diff(range(tmp1)) * 0.3
          }
          else{
            range_var[[variables[1]]][i, j, 1] <- range(tmp1)[1] - diff(range(tmp1)) * 0.3
            range_var[[variables[1]]][i, j, 2] <- range(tmp1)[2] + diff(range(tmp1)) * 0.3

            range_var[[variables[2]]][i, j, 1] <- range(tmp2)[1] - diff(range(tmp2)) * 0.3
            range_var[[variables[2]]][i, j, 2] <- range(tmp2)[2] + diff(range(tmp2)) * 0.3
          }

        }
        # Compute the breaks
        breaks1 <- seq(from = range_var[[variables[1]]][i, j, 1],
                       to = range_var[[variables[1]]][i, j, 2],
                       length.out = nbins1d + 1)

        breaks2 <- seq(from = range_var[[variables[2]]][i, j, 1],
                       to = range_var[[variables[2]]][i, j, 2],
                       length.out = nbins1d + 1)

        # Replace values lower or higher than the range with min or max
        min_range <- range_var[[variables[1]]][i, j, 1]
        max_range <- range_var[[variables[1]]][i, j, 2]
        tmp1[tmp1 < min_range] <- min_range
        tmp1[tmp1 > max_range] <- max_range

        # Replace values lower or higher than the range with min or max
        min_range <- range_var[[variables[2]]][i, j, 1]
        max_range <- range_var[[variables[2]]][i, j, 2]
        tmp2[tmp2 < min_range] <- min_range
        tmp2[tmp2 > max_range] <- max_range

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
  remove(m, v)

  # Return the output as a list of two elements
  output <- list(pdf_matrix, range_var, x_breaks, y_breaks)
  return(output)
}