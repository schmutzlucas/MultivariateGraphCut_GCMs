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

generate_climate_model <- function(n, correlation, temp_mean, temp_sd, precip_mean, precip_sd, humidity_mean, humidity_sd) {
  temp <- rnorm(n, mean = temp_mean, sd = temp_sd)
  precip <- rnorm(n, mean = precip_mean + correlation * (temp - temp_mean), sd = precip_sd)
  humidity <- rnorm(n, mean = humidity_mean + correlation * (temp - temp_mean), sd = humidity_sd)
  joint_dist <- cbind(temp, precip, humidity)
  return(joint_dist)
}



# Generate new model parameters
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

# Generate 4 climate models with 1000 observations each
climate_models <- generate_climate_models(n = 1e4, num_models = 3)

# Print the array
summary(climate_models[,3,1])



