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
NormalizeVariables <- function (data, variables, method) {
  tmp_dim <- dim(data$present)
  tmp_dim[4] <- tmp_dim[4] / length(variables)
  data_nrm <- list()
  data_nrm$future <- data_nrm$present <- array(0, dim = dim(data$present))
  print(tmp_dim)
  for(i in 1:length(variables)) {
    tmp <- Normalize(data$present[, , , i], method)
    data_nrm$present[, , , i] <- array(tmp, dim = tmp_dim)

    tmp <- Normalize(data$future[, , , i], method)
    data_nrm$future[, , , i] <- array(tmp, dim = tmp_dim)
  }
  return(data_nrm)
}