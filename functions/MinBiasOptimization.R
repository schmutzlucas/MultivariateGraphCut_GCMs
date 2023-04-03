#' MinBias optimization method
#'
#' @description
#' Computes the MinBias optimization method that is optimal for the list of variables passed as argument.
#'
#' @param reference A 4D array [lon, lat, 1, n_variables] containing the reference data
#' @param models A 4D array [lon, lat, model, nvariables] containing the model data
#'
#' @return A 2D array with the same size as the input reference data. Each element in the array corresponds to the index of the optimal label for the corresponding location in the reference data.
#'
#' @examples
#' # Load example data
#' data(example_data)
#'
#' # Compute MinBias optimization
#' label_attribution <- MinBiasOptimization(reference = example_data$present[, , 1, ],
#'                                           models = example_data$present[, , 2:4, ])
#'
#' @export
MinBiasOptimization <- function (reference, models) {
  # Number of longitudes
  width <- ncol(reference)
  # Number of latitudes
  height <- nrow(reference)
  # Number of labels used in the GC
  n_labels <- length(models[1, 1, , 1])
  # Number of variables in the reference data
  n_variables <- length(reference[1, 1, 1, ])

  # Initialize bias array
  bias <- array(0, c(height, width, n_labels))
  model_bias <- array(0, c(n_labels))

  # Compute bias for each label and variable
  for (i in 1:n_labels) {
    for (j in 1:n_variables) {
      bias[,,i] <- bias[,,i] + abs(models[,, i, j] - reference[ , , , j])
    }
    model_bias[i] <- sum(bias[,,i])
  }

  # Compute optimal label attribution for each location
  label_attribution <- matrix(0, height, width)
  for(x in 1:height){
    for(y in 1:width){
      label_attribution[x,y] <- which.min(abs(bias[x,y,])) - 1
    }
  }

  # Return optimal label attribution
  return(label_attribution)
}
