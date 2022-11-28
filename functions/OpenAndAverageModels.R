#' @title
#' Opens NETCDF models
#'
#' @description
#'
#'
#' @param model_names
#' @param variables
#'
#' @examples
#'
#' @return
#'
#'
OpenAndAverageModels <- function (model_names, variables,
                                 year_present, year_future) {
  data <- list()
  data_matrix <- list()
  j <- 1
  for(var in variables){
    i <- 1
    for(model_name in model_names){
      # Adjust the filepath and add the necessary suffixes and prefixes
      nc <- nc_open(paste0('data/CMIP5/', var, '/', var, '_', model_name, '.nc'))
      if(j == 1 && i == 1) {
        lat <- ncvar_get(nc, "lat")
        lon <- ncvar_get(nc, "lon")
        data_matrix$present <- data_matrix$future <-
          array(0, c(length(lon), length(lat),
                          length(model_names), length(variables)))
      }

      print(i)
      yyyy <- substr(as.character(nc.get.time.series(nc)), 1, 4)
      iyyyy <- which(yyyy %in% year_present)
      tmp <- apply(
        ncvar_get(nc, var, start = c(1, 1, min(iyyyy)), count = c(-1, -1, length(iyyyy))),
        1:2,
        mean
      )
      data_matrix$present[, , i, j] <- tmp
      data[['present']][[var]][[model_name]] <- tmp

      # future
      iyyyy <- which(yyyy %in% year_future)
      tmp <- apply(
        ncvar_get(nc, var, start = c(1, 1, min(iyyyy)), count = c(-1, -1, length(iyyyy))),
        1:2,
        mean
      )
      data_matrix$future [, , i, j] <- tmp
      data[['future']][[var]][[model_name]] <- tmp
      nc_close(nc)
      i <- i + 1
    }
    j <- j + 1
  }
  remove(j, i)
  # dimnames(data_matrix$present) <- dimnames(data_matrix$future) <-
  #   list(lon = lon, lat = lat, model = paste0(model_names, model_names))

  output <- list(data, data_matrix)
  return(output)
}

