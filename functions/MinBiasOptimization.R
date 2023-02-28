#' MinBias optimization method
#'
#' @description
#' This function produces a labelling using the MinBias optimization method that is optimal for the list of variables passed as argument.
#'
#' @param reference A list of 2D arrays of the reference data, each element of the list being one variable.
#' @param models A list of 3D arrays [lon, lat, model], each element of the list being one variable.
#'
#' @examples
#' # Create reference and model data
#' reference <- list(array(runif(100), dim=c(10,10,1)), array(runif(100), dim=c(10,10,1)))
#' models <- list(array(runif(100), dim=c(10,10,2)), array(runif(100), dim=c(10,10,2)))
#'
#' # Apply the function
#' result <- MinBiasOptimization(reference, models)
#'
#' @return A 2D array with the same size as the input reference data. Each element in the array corresponds to the index of the optimal label for the corresponding location in the reference data.
#'
#' @export
MinBiasOptimization <- function (reference, models) {
  # Number of longitudes
  width <- ncol(reference)
  # Number of latitudes
  height <- nrow(reference)
  # Number of labels used in the GC
  n_labels <- length(models[1,1,,1])
  # Number of variables in the reference data
  n_variables <- length(reference[1, 1, 1, ])

  # Initialize bias array
  bias <- array(0, c(height, width, n_labels))

  # Compute bias for each label and variable
  for (i in 1:n_labels) {
    for (j in 1:n_variables) {
      bias[,,i] <- bias[,,i] + abs(models[,, i, j] - reference[ , , 1, j])
    }
  }

  # Compute optimal label attribution for each location
  label_attribution <- matrix(0, height, width)
  for(x in 1:height){
    for(y in 1:width){
      label_attribution[x,y] <- which.min(abs(bias[x,y,]))-1
    }
  }

  # Return optimal label attribution
  return(label_attribution)
}
