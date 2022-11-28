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
OpenModels <- function(model_names, variables) {
  data <- list()
  for(model_name in model_names) {
    for(var in variables) {
      # Adjust the filepath and add the necessary suffixes and prefixes
      tmp <- nc_open(paste0('data/CMIP5/', var, '/', var, '_',
                            model_name, '.nc'))
      data[[var]][[model_name]] <- ncvar_get(tmp, varid = var)
      rm(tmp)
    }
  }
  return(data)
}