#' @title
#' Normalizes an n dimensions array
#'
#' @description
#' Normalizes an array with the chosen method
#'
#' @param data should be an array
#' @param method allows the user the chose the normalization method.
#' Currently: Standard Score or Min Max
#'
#' @examples
#' Normalize(precipitation, StdSc)
#' Normalize(temperature, MinMax)
#'
#' @return
#' Returns an array of the same size of the one used as argument normalized
#' with the chose method
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