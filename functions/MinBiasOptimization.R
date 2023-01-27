#' @title
#' Min Bias optimization method
#'
#' @description
#' Produces a labelling using the MinBias optimization method
#'
#' @param reference list of 2d arrays of the reference data, each element of the
#'                  list being one variable
#' @param models list of 3d arrays[lon, lat, model], each element of the
#'               list being one variable
#'
#' @examples
#' Normalize(precipitation, StdSc)
#' Normalize(temperature, MinMax)
#'
#' @return
#' Returns an array the same size as the one used as argument normalized
#' with the chosen method
#' 
MinBiasOptimization <- function (reference, models) {
# TODO add other methods
  width <- ncol(reference[[1]])         # Number of longitudes
  height <- nrow(reference[[1]])        # Number of latitudes
  n_labels <- length(models[[1]][1,1,]) # Number of labels used in the GC

  bias <- array(0, c(height, width, nlabs))
  n_variables <- length(reference)

  for(j in 1:n_variables){
    tmp <- models[[j]]
    for(i in 1:n_labels){
      bias[,,i] <- bias[,,i] + (tmp[,,i] - reference[[j]])
    }
  }

  label_attribution <- matrix(0, height, width)
  for(x in 1:height){
    for(y in 1:width){
      label_attribution[x,y] <- which.min(abs(bias[x,y,]))-1
    }
  }
  return(label_attribution)
}