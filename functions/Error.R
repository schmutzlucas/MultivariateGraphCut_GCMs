#' @title
#' Computes the difference between list of arrays
#'
#' @description
#' Error computes the difference between two list of indentically
#' sized array. The lists should be the same length
#'
#' @param ref is a list of length l of arrays of size n*m
#' @param a is a list of length l of arrays of size n*m
#'
#' @return A list of arrays where each cell is the bias between the
#' {a} and {ref}
#'
#' @examples
#' Error(reference_future, labelling_result_future)
#'
Error <- function(ref, a) {
  if (is.list(ref) && is.list(a)){
    error_map <- list(matrix(0, nrow = nrow(as.matrix(ref[[1]])),
                             ncol = ncol(as.matrix(ref[[1]]))))
    l <- length(ref)
    for (i in 1:l){
      error_map[[i]] <- a[[i]] - ref[[i]]
    }
    # TODO implement when ref and a are not list
  } else {
    print('In function Error, ref or a is not a list')
  }
  return(error_map)
}