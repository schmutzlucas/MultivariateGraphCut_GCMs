#' Generate climate models
#'
#' This function generates a given number of synthetic climate models, each consisting of n observations
#' of temperature, precipitation, and humidity variables.
#'
#' @param n The number of observations for each climate model.
#' @param num_models The number of climate models to generate.
#'
#' @return A 3-dimensional array of climate data with dimensions n x 3 x num_models, where n is the
#' number of observations, 3 corresponds to the temperature, precipitation, and humidity variables, and
#' num_models is the number of models.
#'
#' @examples
#' # Generate 3 climate models with n observations each
#' climate_models <- generate_climate_models(n = 10000, num_models = 3)
#'
#' @seealso
#' generate_model_params, generate_climate_model
#'
generate_climate_models <- function(n, num_models) {
  correlation <- 0.5

  # Generate model parameters
  model_params <- generate_model_params(num_models)

  # Create the array to store the models
  climate_models <- array(0, dim = c(n, 3, num_models))

  # Generate the models and store them in the array
  for (i in seq_along(model_params)) {
    model <- generate_climate_model(n = n, correlation = correlation,
                                    temp_mean = model_params[[i]]$temp_mean, temp_sd = model_params[[i]]$temp_sd,
                                    precip_mean = model_params[[i]]$precip_mean, precip_sd = model_params[[i]]$precip_sd,
                                    humidity_mean = model_params[[i]]$humidity_mean, humidity_sd = model_params[[i]]$humidity_sd)
    climate_models[, , i] <- model
  }

  # Return the array of climate models
  return(climate_models)
}


#' Generate synthetic climate data
#'
#' This function generates synthetic climate data consisting of n observations of temperature, precipitation,
#' and humidity variables.
#'
#' @param n The number of observations to generate.
#' @param correlation The correlation between temperature and the other variables.
#' @param temp_mean The mean temperature value.
#' @param temp_sd The standard deviation of the temperature values.
#' @param precip_mean The mean precipitation value.
#' @param precip_sd The standard deviation of the precipitation values.
#' @param humidity_mean The mean humidity value.
#' @param humidity_sd The standard deviation of the humidity values.
#'
#' @return A matrix of size n x 3, where the columns correspond to temperature, precipitation, and humidity.
#'
#' @examples
#' # Generate synthetic climate data with 10000 observations
#' data <- generate_climate_model(n = 10000, correlation = 0.5, temp_mean = 20, temp_sd = 2,
#'                                precip_mean = 50, precip_sd = 10, humidity_mean = 70, humidity_sd = 5)
#'
#' @seealso
#' generate_climate_models, generate_model_params
#'
generate_climate_model <- function(n, correlation, temp_mean, temp_sd, precip_mean, precip_sd, humidity_mean, humidity_sd) {
  temp <- rnorm(n, mean = temp_mean, sd = temp_sd)
  precip <- rnorm(n, mean = precip_mean + correlation * (temp - temp_mean), sd = precip_sd)
  humidity <- rnorm(n, mean = humidity_mean + correlation * (temp - temp_mean), sd = humidity_sd)
  joint_dist <- cbind(temp, precip, humidity)
  return(joint_dist)
}


#' Generate model parameters for synthetic climate data
#'
#' This function generates parameters for multiple synthetic climate models, where each model has its own set of parameters for the mean and standard deviation of the temperature, precipitation, and humidity variables.
#'
#' @param num_models The number of models to generate parameters for.
#'
#' @return A list of length num_models, where each element is a list containing the parameters for one of the synthetic climate models.
#' The list includes the following elements:
#'
#' temp_mean: The mean temperature value.
#'
#' temp_sd: The standard deviation of the temperature values.
#'
#' precip_mean: The mean precipitation value.
#'
#' precip_sd: The standard deviation of the precipitation values.
#'
#' humidity_mean: The mean humidity value.
#'
#' humidity_sd: The standard deviation of the humidity values.
#'
#' @examples
#' # Generate parameters for 5 synthetic climate models
#' params <- generate_model_params(num_models = 5)
#'
#' @seealso
#' generate_climate_models, generate_climate_model
#'
generate_model_params <- function(num_models) {
  temp_mean_range <- c(18, 22)
  temp_sd_range <- c(1.5, 2.5)
  precip_mean_range <- c(50, 60)
  precip_sd_range <- c(10, 15)
  humidity_mean_range <- c(68, 75)
  humidity_sd_range <- c(5, 10)

  temp_means <- seq(from = temp_mean_range[1], to = temp_mean_range[2], length.out = num_models)
  temp_sds <- seq(from = temp_sd_range[1], to = temp_sd_range[2], length.out = num_models)
  precip_means <- seq(from = precip_mean_range[1], to = precip_mean_range[2], length.out = num_models)
  precip_sds <- seq(from = precip_sd_range[1], to = precip_sd_range[2], length.out = num_models)
  humidity_means <- seq(from = humidity_mean_range[1], to = humidity_mean_range[2], length.out = num_models)
  humidity_sds <- seq(from = humidity_sd_range[1], to = humidity_sd_range[2], length.out = num_models)

  model_params <- vector("list", length = num_models)
  for (i in seq_along(model_params)) {
    model_params[[i]] <- list(
      temp_mean = temp_means[i],
      temp_sd = temp_sds[i],
      precip_mean = precip_means[i],
      precip_sd = precip_sds[i],
      humidity_mean = humidity_means[i],
      humidity_sd = humidity_sds[i]
    )
  }
  return(model_params)
}
