#' @title
#' Normalize variables
#'
#' @description
#' This function normalizes an array with the chosen method. It is used to
#' normalize individual variables to make them comparable in the graph cut
#' framework.
#'
#' @param data A list containing the array to be normalized. The array should be
#' of the form [lon, lat, model, var].
#' @param variables A vector of the names of variables present in the array.
#' @param method Allows the user to choose the normalization method.
#' Currently available methods are "StdSc" for standard score and "MinMax" for
#' minimum-maximum normalization.
#'
#' @examples
#' Normalize(data, "tas", "StdSc")
#' Normalize(data, c("tas", "pr"), "MinMax")
#'
#' @return
#' Returns an array of the same size as the one used as argument, but
#' normalized with the chosen method.
#'
NormalizeVariables <- function(data, variables, method) {
  tmp_dim <- dim(data$present)
  data_nrm <- list()
  data_nrm$future <- data_nrm$present <- array(0, dim = tmp_dim)
  for (i in 1:length(variables)) {
    # Normalize all models of one variable together to keep the bias between
    # models.
    tmp <- Normalize(data$present[, , , i], method)
    data_nrm$present[, , , i] <- array(tmp, dim = tmp_dim)

    tmp <- Normalize(data$future[, , , i], method)
    data_nrm$future[, , , i] <- array(tmp, dim = tmp_dim)
  }
  return(data_nrm)
}
