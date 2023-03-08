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


compute_2d_kde <- function(data, var1_idx, var2_idx, nbins, range1, range2) {
  kde_models <- vector("list", length = dim(data)[3])
  for (i in seq_along(kde_models)) {
    kde_models[[i]] <- kde2d(data[, var1_idx, i], data[, var2_idx, i], n = nbins, lims = c(range1, range2))
  }
  return(kde_models)
}

compute_hellinger_dist <- function(kde1, kde2) {
  h_dist <- sqrt(sum((sqrt(kde1$z) - sqrt(kde2$z))^2)) / sqrt(2)
  return(h_dist)
}

compute_all_hellinger_dist <- function(kde_models, kde_ref) {
  num_models <- length(kde_models)
  h_dist_list <- vector("list", length = num_models)
  for (i in seq_along(h_dist_list)) {
    h_dist_list[[i]] <- compute_hellinger_dist(kde_models[[i]], kde_ref)
  }
  return(h_dist_list)
}


# Generate 3 climate models with n observations each
climate_models <- generate_climate_models(n = 1e4, num_models = 10)

reference <- generate_climate_models(n = 1e4, num_models = 1)




# Compute the kde2d objects for all models
pdf_models_list <- compute_kde2d(climate_models, num_models = 3)

# Plot the kde2d object for the third model
ggplot(data.frame(x = kde2d_list[[3]]$x, y = kde2d_list[[3]]$y, z = kde2d_list[[3]]$z), aes(x, y, z = z)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_gradientn(colors = jet.colors(20))


temp_range <- range(reference[,1,])
precip_range <- range(reference[ , 2, ])

# Compute the kde2d objects for all models
tp_kde_models <- compute_2d_kde(climate_models, var1_idx = 1, var2_idx = 2, nbins = 50, range1 = temp_range, range2 = precip_range)
tp_kde_ref <- compute_2d_kde(reference, var1_idx = 1, var2_idx = 2, nbins = 50, range1 = temp_range, range2 = precip_range)

# Compute the Hellinger distance between all models and the reference
h_dist_list <- compute_all_hellinger_dist(tp_kde_models, tp_kde_ref[[1]])

contour(tp_kde_models[[10]])


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
