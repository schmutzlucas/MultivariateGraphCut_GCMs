#' Compute 2D kernel density estimates for given climate data
#'
#' This function computes 2D kernel density estimates for each climate model in the given 3-dimensional array
#' of climate data, using the kde2d() function from the stats package.
#'
#' @param data A 3-dimensional array of climate data with dimensions n x 3 x num_models, where n is the
#' number of observations, 3 corresponds to the temperature, precipitation, and humidity variables, and
#' num_models is the number of models.
#' @param var1_idx The index of the first variable to use for computing the 2D kernel density estimates.
#' This should be an integer between 1 and 3, corresponding to the temperature (1), precipitation (2), or
#' humidity (3) variable.
#' @param var2_idx The index of the second variable to use for computing the 2D kernel density estimates.
#' This should be an integer between 1 and 3, corresponding to the temperature (1), precipitation (2), or
#' humidity (3) variable.
#' @param nbins The number of bins to use for each dimension in the kernel density estimates.
#' @param range1 The range of values to use for the first variable. This should be a numeric vector of length 2.
#' @param range2 The range of values to use for the second variable. This should be a numeric vector of length 2.
#'
#' @return A list of length num_models, where each element is a kde2d() object containing the kernel density
#' estimate for the two variables for one of the climate models.
#'
#' @examples
#' # Generate 3 climate models with n observations each
#' climate_models <- generate_climate_models(n = 10000, num_models = 3)
#'
#' # Compute 2D kernel density estimates for temperature and precipitation variables
#' tp_kde_models <- compute_2d_kde(climate_models, var1_idx = 1, var2_idx = 2, nbins = 50,
#' range1 = range(climate_models[, 1, ]), range2 = range(climate_models[, 2, ]))
#'
#' # Compute 2D kernel density estimates for temperature and humidity variables
#' th_kde_models <- compute_2d_kde(climate_models, var1_idx = 1, var2_idx = 3, nbins = 50,
#' range1 = range(climate_models[, 1, ]), range2 = range(climate_models[, 3, ]))
#'
#' @seealso
#' kde2d()
#'
#' @importFrom stats kde2d
compute_2d_kde <- function(data, var1_idx, var2_idx, nbins, range1, range2) {
  kde_models <- vector("list", length = dim(data)[3])
  for (i in seq_along(kde_models)) {
    kde_models[[i]] <- kde2d(data[, var1_idx, i], data[, var2_idx, i], n = nbins, lims = c(range1, range2))
  }
  return(kde_models)
}


#' Compute Hellinger distance between two 2D kernel density estimates
#'
#' Given two 2D kernel density estimates computed using the kde2d() function, this function
#' computes the Hellinger distance between them.
#'
#' @param kde1 A kde2d() object containing the kernel density estimate for the first variable.
#' @param kde2 A kde2d() object containing the kernel density estimate for the second variable.
#'
#' @return A scalar value representing the Hellinger distance between the two kernel density estimates.
#'
#' @examples
#' # Generate 3 climate models with n observations each
#' climate_models <- generate_climate_models(n = 10000, num_models = 3)
#'
#' # Compute 2D kernel density estimates for temperature and precipitation variables
#' tp_kde_models <- compute_2d_kde(climate_models, var1_idx = 1, var2_idx = 2, nbins = 50,
#'                                 range1 = range(climate_models[, 1, ]), range2 = range(climate_models[, 2, ]))
#'
#' # Compute the Hellinger distance between the first and second climate models
#' h_dist <- compute_hellinger_dist(tp_kde_models[[1]], tp_kde_models[[2]])
#'
#' @seealso
#' kde2d
#'
#' @importFrom stats kde2d
compute_hellinger_dist <- function(kde1, kde2) {
  h_dist <- sqrt(sum((sqrt(kde1$z) - sqrt(kde2$z))^2)) / sqrt(2)
  return(h_dist)
}

#' Compute Hellinger distances between multiple 2D kernel density estimates and a reference model
#'
#' Given a list of 2D kernel density estimates and a reference 2D kernel density estimate, this function
#' computes the Hellinger distance between each of the estimates and the reference.
#'
#' @param kde_models A list of kde2d() objects containing the kernel density estimates for each of
#' the climate models.
#' @param kde_ref A kde2d() object containing the kernel density estimate for the reference model.
#'
#' @return A list of scalar values representing the Hellinger distance between each of the kernel density
#' estimates and the reference model.
#'
#' @examples
#' # Generate 10 climate models with n observations each
#' climate_models <- generate_climate_models(n = 10000, num_models = 10)
#'
#' # Generate a reference climate model with n observations
#' reference <- generate_climate_models(n = 10000, num_models = 1)
#'
#' # Compute 2D kernel density estimates for temperature and precipitation variables for all models and the reference
#' tp_kde_models <- compute_2d_kde(climate_models, var1_idx = 1, var2_idx = 2, nbins = 50,
#'                                 range1 = range(reference[, 1, ]), range2 = range(reference[, 2, ]))
#' tp_kde_ref <- compute_2d_kde(reference, var1_idx = 1, var2_idx = 2, nbins = 50,
#'                              range1 = range(reference[, 1, ]), range2 = range(reference[, 2, ]))
#'
#' # Compute the Hellinger distance between all models and the reference
#' h_dist_list <- compute_all_hellinger_dist(tp_kde_models, tp_kde_ref[[1]])
#'
#' @seealso
#' kde2d
#'
#' @importFrom stats kde2d
compute_all_hellinger_dist <- function(kde_models, kde_ref) {
  num_models <- length(kde_models)
  h_dist_list <- vector("list", length = num_models)
  for (i in seq_along(h_dist_list)) {
    h_dist_list[[i]] <- compute_hellinger_dist(kde_models[[i]], kde_ref)
  }
  return(h_dist_list)
}