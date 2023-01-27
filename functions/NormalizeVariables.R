#' @title
#' Normalize variables
#'
#' @description
#' Normalizes an array with the chosen method.
#' Used to normalized individual variables to make them comparable in the graph
#' cut framework
#'
#' @param data list containing the array to be normalized, e.g. list(present, future)
#'             the array sould be [lon, lat, model, var]
#' @param variables vector of the names of variables present in the array
#' @param method allows the user the chose the normalization method.
#' Currently: Standard Score or Min Max
#'
#' @examples
#' Normalize(precipitation, StdSc)
#' Normalize(temperature, MinMax)
#'
#' @return
#' Returns an array of the same size of the one used as argument normalized
#' with the chosen method
#'
NormalizeVariables <- function (data, variables, method) {
  tmp_dim <- dim(data[1])
  tmp_dim[4] <- tmp_dim[4] / length(variables)
  data_nrm <- list()
  data_nrm$future <- data_nrm$present <- array(0, dim = dim(data$present))
  print(tmp_dim)
  for(i in 1:length(variables)) {
    #TODO Should we normalize data model by model, or all models together?
    # Now we normalize all models of one variables togeter => the bias between
    # models is kept
    tmp <- Normalize(data$present[, , , i], method)
    data_nrm$present[, , , i] <- array(tmp, dim = tmp_dim)

    tmp <- Normalize(data$future[, , , i], method)
    data_nrm$future[, , , i] <- array(tmp, dim = tmp_dim)
  }
  return(data_nrm)
}