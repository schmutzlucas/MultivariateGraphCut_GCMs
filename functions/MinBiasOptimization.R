#' @title
#' Min Bias optimization method
#'
#' @description
#' Produces a labelling using the MinBias optimization method
#'
#' @param referenceerence list of 2d arrays of the referenceerence data, each entry being one modelsiable
#' @param models list of 3d arrays
#'
#' @examples
#' Normalize(precipitation, StdSc)
#' Normalize(temperature, MinMax)
#'
#' @return
#' Returns an array of the same size of the one used as argument normalized
#' with the chose method
#' 
MinBiasOptimization <- function (reference, models) {
# TODO
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