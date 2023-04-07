#' @title
#' Opens NETCDF models and averages variables over time
#' @description
#' This function reads and extracts variables from NETCDF files in a specific
#' directory and calculates the average values over a given time period for
#' present and future. The result is returned as a list of two elements: a list
#' containing the extracted data and a list containing the average values for
#' each variable and model.
#'
#' @param model_names a character vector of model names, e.g., c("model1", "model2")
#' @param variables a character vector of variable names, e.g., c("temp", "precip")
#' @param year_present a numeric vector specifying the years for present, e.g., c(1981, 2010)
#' @param year_future a numeric vector specifying the years for future, e.g., c(2071, 2100)
#'
#' @return a list containing two elements:
#' - data: a list containing the extracted data for each variable and model
#' - data_matrix: a list containing the average values for each variable and model
#'
#' @examples
#' # Extract temperature and precipitation data for present and future time periods from two models
#' data <- OpenAndAverageModels(c("model1", "model2"), c("temp", "precip"), c(1981, 2010), c(2071, 2100))
#'
#' # Print the extracted data for the first model and the variable "temp"
#' print(data$data$present$temp$model1)
#'
#' # Print the average values for the variable "precip" and all models
#' print(data$data_matrix$present[, , , "precip"])
#'
#' @import ncdf4
#'
OpenAndAverageCMIP6 <- function (model_names, variables,
                                 year_present, year_future, period) {

  # Initialize data structures
  data <- list()
  data_matrix <- list()

  # Loop through variables and models
  j <- 1
  for(var in variables){
    i <- 1
    for(model_name in model_names){
      dir_path <- paste0('data/CMIP6/', model_name, '/', var, '/')
      # Create the pattern
      pattern <- glob2rx(paste0(var, "_", model_name, "_", period, "*.nc"))

      # Get the filepath
      file_name <- list.files(path = dir_path,
                             pattern = pattern)
      file_path <- paste0(dir_path, file_name)
      print(file_path)

      # Check that there is only one matching file
      if (length(file_path) == 1) {
        nc <- nc_open(file_path)
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
      } else {
        # Handle the case where there are multiple or no matching files
        print("Error: Found multiple or no matching files")
      }

      

      # Update counter for models
      i <- i + 1
    }

    # Update counter for variables
    j <- j + 1
  }

  # Remove the counters
  remove(j, i)

  # Set the dimnames of the matrices
  dimnames(data_matrix$present) <- dimnames(data_matrix$future) <-
    list(lon = lon, lat = lat, model = paste0(model_names, model_names))

  # Return the output as a list of two elements
  output <- list(data, data_matrix)
  return(output)
}

