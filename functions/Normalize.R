#' @title
#' Normalize array
#'
#' @description
#' Normalizes an array with the chosen method.
#'
#' @param data array to be normalized
#' @param method normalization method to be used. Currently available methods are 'StdSc'
#' and 'min_max'
#'
#' @examples
#' Normalize(precipitation, method = 'StdSc')
#' Normalize(temperature, method = 'min_max')
#'
#' @return
#' Returns the normalized array of the same size as the input array.
#'
Normalize <- function(data, method){

  if (method == 'StdSc') {
    data <- (data - mean(data)) / (sd(data))
  }
  else if (method == 'min_max') {
    data <- (data - min(data)) / (max(data) - min(data))
  }
  else {
    print('method argument is invalid')
  }
  return(data)
}
