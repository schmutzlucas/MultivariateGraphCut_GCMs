#' Opens NETCDF models
#'
#' This function opens and reads in NETCDF model files for a given set of models and variables.
#' It returns the resulting data in a list structure, where each variable has a nested list of model data.
#'
#' @param model_names A character vector representing the names of the models to open and read.
#' @param variables A character vector representing the names of the variables to read from the models.
#'
#' @examples
#' # Load example data
#' data("example_model_names")
#' data("example_variables")
#'
#' # Open models and read data
#' models_data <- OpenModels(
#'   model_names = example_model_names,
#'   variables = example_variables
#' )
#'
#' # Print data dimensions
#' for(var in example_variables) {
#'   for(model_name in example_model_names) {
#'     cat(paste(var, model_name, dim(models_data[[var]][[model_name]])), "\n")
#'   }
#' }
#'
#' @return A list containing the data for the specified models and variables.
#'
#' @references
#' https://www.unidata.ucar.edu/software/netcdf/
#'
#' @import ncdf4
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
