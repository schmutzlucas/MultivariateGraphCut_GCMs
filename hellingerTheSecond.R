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

compute_joint_dist <- function(model, var1_idx, var2_idx, nbins) {
  # Extract the data for the two variables
  X_data <- model[, var1_idx]
  Y_data <- model[, var2_idx]

  # Compute the joint histogram counts for the two variables
  XY_hist <- hist2d(X_data, Y_data, nbins = nbins)

  # Compute the joint probabilities for the two variables
  XY_prob <- XY_hist$counts / sum(XY_hist$counts)

  # Return the joint probabilities
  return(XY_prob)
}

compute_kde2d <- function(climate_models, num_models, n_bins = 50) {
  kde2d_list <- vector("list", length = num_models)

  temp_range <- range(climate_models[, 1, ])
  precip_range <- range(climate_models[, 2, ])

  for (i in seq_along(kde2d_list)) {
    kde2d_list[[i]] <- kde2d(climate_models[, 1, i], climate_models[, 2, i],
                             n = n_bins, lims = c(temp_range, precip_range))
  }

  return(kde2d_list)
}


# Generate 3 climate models with n observations each
climate_models <- generate_climate_models(n = 1e4, num_models = 3)

# Generate the reference model
reference <- generate_climate_models(n = 1e4, num_models = 1)[, , 1]

# Compute the range of temperature and precipitation across both models
temp_range <- range(climate_models[, 1, ])
precip_range <- range(climate_models[, 2, ])

# Compute the 2D kernel density estimate of temperature and precipitation for the second model
tp_dens_model2 <- kde2d(climate_models[, 1, 2], climate_models[, 2, 2], n = 50,
                        lims = c(temp_range, precip_range))

# Compute the 2D kernel density estimate of temperature and precipitation for the third model
tp_dens_model3 <- kde2d(climate_models[, 1, 3], climate_models[, 2, 3], n = 50,
                       lims = c(temp_range, precip_range))






h_dist <- sqrt(sum((sqrt(tp_dens_model3$z) - sqrt(tp_dens_model2$z))^2)) / sqrt(2)











contour(tp_dens_model3)


# Convert the density matrix to long format for ggplot
# Convert the density matrix to long format for ggplot
tp_dens_model3_long <- reshape2::melt(tp_dens_model3$z)

# Create a ggplot graph of the density estimates
ggplot(tp_dens_model3_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_viridis() +
  labs(x = "Temperature", y = "Precipitation", fill = "Density") +
  ggtitle("Density Estimation for Climate Model 3") +
  theme_bw()

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

# Create a ggplot graph of the density estimates with jet color scale
ggplot(tp_dens_model3_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = jet.colors(20)) +
  labs(x = "Temperature", y = "Precipitation", fill = "Density") +
  ggtitle("Density Estimation for Climate Model 3") +
  theme_bw()
