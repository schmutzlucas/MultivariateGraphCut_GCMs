# Install and load necessary libraries
list_of_packages <- read.table("package_list.txt", sep="\n")$V1
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages))
  install.packages(new.packages, repos = "https://cloud.r-project.org")

library(devtools)
lapply(list_of_packages, library, character.only = TRUE)
# install_github("schmutzlucas/gcoWrapR")

# Loading local functions
source_code_dir <- 'functions/'  # The directory where all functions are saved.
file_paths <- list.files(source_code_dir, full.names = TRUE)
for(path in file_paths){
  source(path)
}


# # Setting global variables
# lon <- -10:10
# lat <- -10:10
# lon_size <- length(lon)
# lat_size <- length(lat)
#
# # Temporal ranges
# year_present <<- 1950:1975
# year_future <<- 1998:2023



# Setting global variables
lon <- -180:179
lat <- -90:90
lon_size <- length(lon)
lat_size <- length(lat)

# Temporal ranges
year_present <<- 1950:1975
year_future <<- 1998:2023

# Data directory
data_dir <<- 'data/CMIP6_merged_all/'

# List of variables used
variables <- c('pr', 'tas', 'psl')

# Bins for the PDFs (nbins1d is the number of bins per variable)
nbins1d <<- 8
nbins <<- nbins1d^(length(variables))
# (The joint PDF will have nbins1d^n_vars bins)

# Obtain the list of models from a file
model_names <- read.table('model_names_pr_tas_psl.txt')
model_names <- as.list(model_names[['V1']])
ref_index <<- 1

# Separate the reference model from the rest.
reference_name <<- model_names[[ref_index]]
model_names <<- model_names[-ref_index]

# Custom function to format time into human-readable format
format_time <- function(time_seconds) {
  hours <- floor(time_seconds / 3600)
  minutes <- floor((time_seconds %% 3600) / 60)
  seconds <- round(time_seconds %% 60, 2)
  if (hours > 0) {
    return(paste(hours, "hours", minutes, "minutes", seconds, "seconds"))
  } else if (minutes > 0) {
    return(paste(minutes, "minutes", seconds, "seconds"))
  } else {
    return(paste(seconds, "seconds"))
  }
}

# # Time the execution of the new bias-corrected function.
# time_optimized <- system.time({
#   results <- compute_nd_pdf_bias_corrected(
#     variables, reference_name, model_names, data_dir,
#     year_present, year_future, lon, lat, nbins1d,
#     range_var = aperm(abind(range_var_final$ranges, along = 4), c(1, 2, 4, 3)), verbose = TRUE
#   )
#
# })
# cat("Time taken for compute_nd_pdf_bias_corrected: ",
#     format_time(time_optimized["elapsed"]), "\n")


# Time the execution of the new bias-corrected function.
time_optimized <- system.time({
  results <- compute_nd_pdf_bias_corrected_2(
    variables,
    reference_name,
    model_names,
    data_dir,
    year_present,
    year_future,
    lon,
    lat,
    nbins1d,
    workers = 4,    # Adjust the number of workers as needed
    buffer = 0.10,
    verbose = TRUE
  )
})
cat("Time taken for compute_nd_pdf_bias_corrected: ",
    format_time(time_optimized["elapsed"]), "\n")



# Extract PDFs from the unified function results
pdf_ref_present <- results$pdf_ref$present
pdf_models_present <- results$pdf_models$present
pdf_ref_future <- results$pdf_ref$future
pdf_models_future <- results$pdf_models$future


# Get the current date and time
current_time <- Sys.time()

# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")

# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "_my_workspace_ERA5_bias_corrected_new_22models_90-90.RData")

# Save the workspace using the generated filename
save.image(file = filename, compress = FALSE)


# Plotting the 3D PDF of one gridpoint
{
  library(plotly)
  library(abind)

  # Example indices (adjust as needed)
  lon_index <- 90 + 1  # selected longitude index
  lat_index <- 90 + 45 # selected latitude index
  model_idx <- 2  # select one model (from the pdf_models_future array)
  nbins <- 8      # number of bins per variable

  # Extract the PDF vector for the chosen grid point and model from the results.
  pdf_vector <- results$pdf_reference$future[lon_index, lat_index, ]

  # Reshape the 1D PDF vector into a 3D array.
  pdf_3d <- array(pdf_vector, dim = c(nbins, nbins, nbins))

  # Adjust range_var to extract the correct reference range for the chosen grid point
  range_mat <- results$ref_range_future[lon_index, lat_index, , ]  # dimensions: [3, 2]
  print(range_mat)


  # Compute bin edges and centers for each variable.
  centers <- list()
  for (v in 1:3) {
    bin_edges <- seq(range_mat[v, 1], range_mat[v, 2], length.out = nbins + 1)
    centers[[v]] <- (bin_edges[-1] + bin_edges[-length(bin_edges)]) / 2
  }

  # Create a grid of bin centers.
  grid <- expand.grid(pr = centers[[1]], tas = centers[[2]], psl = centers[[3]])

  # Flatten the 3D PDF into a vector.
  pdf_flat <- as.vector(pdf_3d)
  # Normalize PDF values for marker size.
  normalized_pdf <- pdf_flat / max(pdf_flat, na.rm = TRUE)

  # Plot the 3D scatter plot
  fig <- plot_ly(
    data = grid,
    x = ~pr,
    y = ~tas,
    z = ~psl,
    type = "scatter3d",
    mode = "markers",
    marker = list(
      size = ~normalized_pdf * 75,  # Adjust scaling factor as needed
      color = ~pdf_flat,
      colorscale = "Viridis",
      showscale = TRUE
    ),
    text = ~paste("PDF Value:", round(pdf_flat, 4))
  ) %>% layout(
    scene = list(
      xaxis = list(title = "pr"),
      yaxis = list(title = "tas"),
      zaxis = list(title = "psl")
    ),
    title = paste("3D PDF for Model", model_idx, "Grid Point (Lon:", lon_index, ", Lat:", lat_index, ")")
  )

  fig
}

# Plotting the 3d pdf of one gridpoint
{
  if (exists("fig")) rm(fig)
  # Example indices (adjust as needed)
  lon_index <- 90 + 1  # selected longitude index
  lat_index <- 90 + 45 # selected latitude index
  model_idx <- 1 # select one model (from the pdf_models array)
  nbins <- 8       # number of bins per variable

  # Extract the PDF vector for the chosen grid point and model.
  # Here pdf_models is from results$pdf_models$present and has dimensions:
  # [lon, lat, nbins^n_vars, num_models].
  pdf_vector <- MMM_future[lon_index, lat_index, ]

  # Reshape the 1D PDF vector into a 3D array.
  pdf_3d <- array(pdf_vector, dim = c(nbins, nbins, nbins))

  # Extract the reference range for that grid point.
  # We assume results$ref_range_present has dimensions: [lon, lat, n_vars, 2]
  range_mat <- results$ref_range_future[lon_index, lat_index, , ]  # dimensions: [3, 2]
  print(range_mat)

  # Compute bin edges and centers for each variable.
  centers <- list()
  for(v in 1:3) {
    bin_edges <- seq(range_mat[v, 1], range_mat[v, 2], length.out = nbins + 1)
    centers[[v]] <- (bin_edges[-1] + bin_edges[-length(bin_edges)])/2
  }

  # Create a grid of bin centers.
  grid <- expand.grid(x = centers[[1]], y = centers[[2]], z = centers[[3]])

  # Flatten the 3D PDF into a vector.
  pdf_flat <- as.vector(pdf_3d)
  # Normalize PDF values for marker size.
  normalized_pdf <- pdf_flat / max(pdf_flat, na.rm = TRUE)

  library(plotly)
  fig <- plot_ly(
    data = grid,
    x = ~x,
    y = ~y,
    z = ~z,
    type = "scatter3d",
    mode = "markers",
    marker = list(
      size = ~normalized_pdf * 75,
      color = ~pdf_flat,
      colorscale = "Viridis",
      showscale = TRUE
    ),
    text = ~paste("PDF Value:", round(pdf_flat, 4))
  ) %>%
    layout(
      scene = list(
        xaxis = list(title = "pr"),
        yaxis = list(title = "tas"),
        zaxis = list(title = "psl"),
        aspectmode = "cube"  # <-- This ensures all axes are treated equally
      ),
      title = paste("3D PDF for Model", model_idx, "Grid Point (Lon:", lon_index, ", Lat:", lat_index, ")")
    )

  fig


}


# Plotting the 3d pdf of one gridpoint
{
  if (exists("fig")) rm(fig)
  # Example indices (adjust as needed)
  lon_index <- 90 + 1  # selected longitude index
  lat_index <- 90 + 45 # selected latitude index
  model_idx <- 1 # select one model (from the pdf_models array)
  nbins <- 8       # number of bins per variable

  # Extract the PDF vector for the chosen grid point and model.
  # Here pdf_models is from results$pdf_models$present and has dimensions:
  # [lon, lat, nbins^n_vars, num_models].
  pdf_vector <- pdf_models_future[lon_index, lat_index, , model_idx]

  # Reshape the 1D PDF vector into a 3D array.
  pdf_3d <- array(pdf_vector, dim = c(nbins, nbins, nbins))

  # Extract the reference range for that grid point.
  # We assume results$ref_range_present has dimensions: [lon, lat, n_vars, 2]
  range_mat <- results$ref_range_future[lon_index, lat_index, , ]  # dimensions: [3, 2]
  print(range_mat)

  # Compute bin edges and centers for each variable.
  centers <- list()
  for(v in 1:3) {
    bin_edges <- seq(range_mat[v, 1], range_mat[v, 2], length.out = nbins + 1)
    centers[[v]] <- (bin_edges[-1] + bin_edges[-length(bin_edges)])/2
  }

  # Create a grid of bin centers.
  grid <- expand.grid(x = centers[[1]], y = centers[[2]], z = centers[[3]])

  # Flatten the 3D PDF into a vector.
  pdf_flat <- as.vector(pdf_3d)
  # Normalize PDF values for marker size.
  normalized_pdf <- pdf_flat / max(pdf_flat, na.rm = TRUE)

  library(plotly)
  fig <- plot_ly(
    data = grid,
    x = ~x,
    y = ~y,
    z = ~z,
    type = "scatter3d",
    mode = "markers",
    marker = list(
      size = ~normalized_pdf * 75,
      color = ~pdf_flat,
      colorscale = "Viridis",
      showscale = TRUE
    ),
    text = ~paste("PDF Value:", round(pdf_flat, 4))
  ) %>%
    layout(
      scene = list(
        xaxis = list(title = "pr"),
        yaxis = list(title = "tas"),
        zaxis = list(title = "psl"),
        aspectmode = "cube"  # <-- This ensures all axes are treated equally
      ),
      title = paste("3D PDF for Model", model_idx, "Grid Point (Lon:", lon_index, ", Lat:", lat_index, ")")
    )

  fig


}

# --- Begin Complete Hellinger Distance Computation using compute_partial_hdist ---

# Create a nested list of selected indices for every grid point (using all bins)
all_bins <- 1:nbins
selected_indices_all <- vector("list", length(lon))
for(i in seq_along(lon)) {
  selected_indices_all[[i]] <- vector("list", length(lat))
  for(j in seq_along(lat)) {
    selected_indices_all[[i]][[j]] <- all_bins
  }
}

# Compute complete Hellinger distances for present and future using compute_partial_hdist:
h_dist_present <- compute_partial_hdist(pdf_ref_present, pdf_models_present, selected_indices_all)
h_dist_future  <- compute_partial_hdist(pdf_ref_future, pdf_models_future, selected_indices_all)

# Plot histograms for debugging
hist(h_dist_present, main = "Complete Hellinger Distance (Present)")
hist(h_dist_future, main = "Complete Hellinger Distance (Future)")

# --- End Complete Hellinger Distance Computation ---

smooth_cost <- 0.6
tryCatch({
  GC_result06 <- GraphCutHellinger_nD_lat(
    pdf_models_future = pdf_models_future,
    h_dist = h_dist_present,
    weight_data = 1,               # Fixed data weight
    weight_smooth = smooth_cost,   # Varying smooth cost
    nBins = nbins1d^3,
    lat = lat,
    seed = 1,
    verbose = TRUE,
    rebuild = TRUE
  )
}, error = function(e) {
  cat("Error encountered with smooth cost =", smooth_cost, ": ", e$message, "\n")
})
gc()

smooth_cost <- 0.1
tryCatch({
  GC_result01 <- GraphCutHellinger_nD_lat(
    pdf_models_future = pdf_models_future,
    h_dist = h_dist_present,
    weight_data = 1,               # Fixed data weight
    weight_smooth = smooth_cost,   # Varying smooth cost
    nBins = nbins1d^3,
    lat = lat,
    seed = 1,
    verbose = TRUE,
    rebuild = TRUE
  )
}, error = function(e) {
  cat("Error encountered with smooth cost =", smooth_cost, ": ", e$message, "\n")
})
gc()

smooth_cost <- 0.2
tryCatch({
  GC_result02 <- GraphCutHellinger_nD_lat(
    pdf_models_future = pdf_models_future,
    h_dist = h_dist_present,
    weight_data = 1,               # Fixed data weight
    weight_smooth = smooth_cost,   # Varying smooth cost
    nBins = nbins1d^3,
    lat = lat,
    seed = 1,
    verbose = TRUE,
    rebuild = TRUE
  )
}, error = function(e) {
  cat("Error encountered with smooth cost =", smooth_cost, ": ", e$message, "\n")
})
gc()

# Map of labels
{
  # Extract the label attribution for the current smooth cost
  GC_labels <- GC_result01$label_attribution

  # Convert the label matrix to a data frame for plotting
  label_df <- reshape2::melt(GC_labels, varnames = c("lon_idx", "lat_idx"), value.name = "label_attribution")

  # Explicitly assign correct longitude and latitude values
  label_df$lon <- lon[label_df$lon_idx]  # Map longitude indices to values
  label_df$lat <- lat[label_df$lat_idx]    # Map latitude indices to values

  # Force label_attribution to be a factor with levels in the exact order of model_names.
  label_df$label_attribution <- factor(label_df$label_attribution,
                                       levels = seq_along(model_names),
                                       labels = model_names)

  # Create a named color palette: each model name is explicitly mapped to its color.
  color_palette <- c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
    "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
    "#393b79", "#5254a3", "#6b6ecf"
  )
  # Name the palette vector with model_names (in the same order)
  names(color_palette) <- model_names

  p6 <- ggplot() +
    geom_tile(data = label_df, aes(x = lon, y = lat, fill = label_attribution)) +
    scale_fill_manual(
      values = color_palette,
      na.value = "white",
      guide = guide_legend(title = "Model Names", ncol = 1)
    ) +
    ggtitle(paste("Label GC Hellinger - Lambda:", 0.1)) +
    borders("world", colour = 'black', size = 0.12) +
    theme_bw() +
    theme(
      legend.position = 'right',
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      legend.key.size = unit(0.5, 'cm'),
      legend.key.height = unit(0.5, 'cm'),
      legend.key.width = unit(0.5, 'cm'),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      plot.title = element_text(size = 16),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    ) +
    xlab('Longitude') +
    ylab('Latitude') +
    easy_center_title()

  p6

  # Generate file name based on the smooth cost
  name <- paste0("figure/GC_labelling_smooth01_BC_22model")

  # Save the plot as both PDF and PNG
  ggsave(paste0(name, ".pdf"), plot = p6, width = 20, height = 15, units = "cm", dpi = 300)
  ggsave(paste0(name, ".png"), plot = p6, width = 20, height = 15, units = "cm", dpi = 300)
}


num_bins <- 20
x_limits <- range(h_dist_present, finite = TRUE)  # Adjust based on your data
y_limits <- c(0, 30000)

# Loop through models and plot histograms
for (v in seq_along(model_names)) {
  hist(
    h_dist_present[,,v],
    breaks = seq(x_limits[1], x_limits[2], length.out = num_bins + 1),
    xlim = x_limits,
    ylim = y_limits,
    main = paste("Histogram for", model_names[v]),
    xlab = "Value",
    col = "lightblue",
    border = "black"
  )
}
# Extract the data
nan_map <- is.na(results$out_of_range_counts$present)

# Create a plot
image(nan_map, col = c("white", "black"), axes = FALSE, main = "Map of NaN Values")

# Add a legend
legend("topright", legend = c("Valid", "NaN"), fill = c("white", "black"), border = "black")



# H Dist GC
# Initialize a lon x lat matrix for each smooth cost
GC01_hdist_present <- matrix(NA, nrow = length(lon), ncol = length(lat))
GC01_hdist_future <- matrix(NA, nrow = length(lon), ncol = length(lat))

GC06_hdist_present <- matrix(NA, nrow = length(lon), ncol = length(lat))
GC06_hdist_future <- matrix(NA, nrow = length(lon), ncol = length(lat))

for(l in 1:(length(model_names))){  # Ensure that indexing aligns with model names
  islabel <- which(GC_result01$label_attribution == l)
  GC01_hdist_present[islabel] <- h_dist_present[,,l][islabel]
  GC01_hdist_future[islabel] <- h_dist_future[,,l][islabel]

  islabel <- which(GC_result06$label_attribution == l)
  GC06_hdist_present[islabel] <- h_dist_present[,,l][islabel]
  GC06_hdist_future[islabel] <- h_dist_future[,,l][islabel]
}



# MMM

# Compute Multi-Model Mean for Present
MMM_present <- apply(pdf_models_present, c(1, 2, 3), mean)

# Compute Multi-Model Mean for Future
MMM_future <- apply(pdf_models_future, c(1, 2, 3), mean)


# Initialize arrays to store the Hellinger distance for MMM
MMM_hdist <- array(NA, dim = c(length(lon), length(lat)))
MMM_hdist_future <- array(NA, dim = c(length(lon), length(lat)))

# Compute Hellinger distance for the Multi-Model Mean
for (i in seq_along(lon)) {
  for (j in seq_along(lat)) {
    # Compute Hellinger distance for present
    MMM_hdist_present[i, j] <- sqrt(sum((sqrt(MMM_present[i, j, ]) - sqrt(pdf_ref_present[i, j, ]))^2)) / sqrt(2)

    # Compute Hellinger distance for future
    MMM_hdist_future[i, j] <- sqrt(sum((sqrt(MMM_future[i, j, ]) - sqrt(pdf_ref_future[i, j, ]))^2)) / sqrt(2)
  }
}
gc()
# Get the current date and time
current_time <- Sys.time()

# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")

# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "_my_workspace_ERA5_bias_corrected_22models_H_dist.RData")

# Save the workspace using the generated filename
save.image(file = filename, compress = FALSE)


# Map H Dist future
test_df <- melt(GC01_hdist_future, c("lon", "lat"), value.name = "H_dist")

p6 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon-180, y=lat-90, fill=H_dist))+
  labs(subtitle = 'Smooth = 0.1, Projection period : 1998 - 2023')+
  ggtitle(paste0('GC Bias corrected', ': Average H = ', round(mean(GC01_hdist_future), 2)))+
  scale_fill_gradient(low = "white", high = "#015a8c", limits = c(0.1, 0.70), oob = scales::squish)+
  borders("world", colour = 'black', lwd = 0.12) +
  scale_x_continuous(, expand = c(0, 0)) +
  scale_y_continuous(, expand = c(0, 0))+
  theme(legend.position = 'bottom')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())+
  xlab('Longitude')+
  ylab('Latitude') +
  labs(fill='Hellinger \nDistance')+
  theme_bw()+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=16), #change legend title font sizen
        legend.text = element_text(size=12))+ #change legend text font size
  theme(plot.title = element_text(size=24),
        plot.subtitle = element_text(size = 20,hjust=0.5),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),)+
  easy_center_title()
p6

# Generate file name based on the smooth cost
name <- paste0("figure/GC_H_dist_smooth01_BC_projection")

# Save the plot as both PDF and PNG
ggsave(paste0(name, ".pdf"), plot = p6, width = 20, height = 15, units = "cm", dpi = 300)
ggsave(paste0(name, ".png"), plot = p6, width = 20, height = 15, units = "cm", dpi = 300)



# Map H Dist future


test_df <- melt(GC06_hdist_future, c("lon", "lat"), value.name = "H_dist")

p6 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon-180, y=lat-90, fill=H_dist))+
  labs(subtitle = 'Smooth = 0.6, Projection period : 1998 - 2023')+
  ggtitle(paste0('GC Bias corrected', ': Average H = ', round(mean(GC06_hdist_future), 2)))+
  scale_fill_gradient(low = "white", high = "#015a8c", limits = c(0.1, 0.70), oob = scales::squish)+
  borders("world", colour = 'black', lwd = 0.12) +
  scale_x_continuous(, expand = c(0, 0)) +
  scale_y_continuous(, expand = c(0, 0))+
  theme(legend.position = 'bottom')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())+
  xlab('Longitude')+
  ylab('Latitude') +
  labs(fill='Hellinger \nDistance')+
  theme_bw()+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=16), #change legend title font sizen
        legend.text = element_text(size=12))+ #change legend text font size
  theme(plot.title = element_text(size=24),
        plot.subtitle = element_text(size = 20,hjust=0.5),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),)+
  easy_center_title()
p6

# Generate file name based on the smooth cost
name <- paste0("figure/GC_H_dist_smooth06_BC_projection")

# Save the plot as both PDF and PNG
ggsave(paste0(name, ".pdf"), plot = p6, width = 20, height = 15, units = "cm", dpi = 300)
ggsave(paste0(name, ".png"), plot = p6, width = 20, height = 15, units = "cm", dpi = 300)


# Map H Dist present


test_df <- melt(GC01_hdist_present, c("lon", "lat"), value.name = "H_dist")

p6 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon-180, y=lat-90, fill=H_dist))+
  labs(subtitle = 'Smooth = 0.1, Calibration period : 1950 - 1975')+
  ggtitle(paste0('GC BC', ': Average H = ', round(mean(GC01_hdist_present), 2)))+
  scale_fill_gradient(low = "white", high = "#015a8c", limits = c(0.1, 0.70), oob = scales::squish)+
  borders("world", colour = 'black', lwd = 0.12) +
  scale_x_continuous(, expand = c(0, 0)) +
  scale_y_continuous(, expand = c(0, 0))+
  theme(legend.position = 'bottom')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())+
  xlab('Longitude')+
  ylab('Latitude') +
  labs(fill='Hellinger \nDistance')+
  theme_bw()+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=16), #change legend title font sizen
        legend.text = element_text(size=12))+ #change legend text font size
  theme(plot.title = element_text(size=24),
        plot.subtitle = element_text(size = 20,hjust=0.5),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),)+
  easy_center_title()
p6

# Generate file name based on the smooth cost
name <- paste0("figure/GC_H_dist_smooth01_BC_calibration")

# Save the plot as both PDF and PNG
ggsave(paste0(name, ".pdf"), plot = p6, width = 20, height = 15, units = "cm", dpi = 300)
ggsave(paste0(name, ".png"), plot = p6, width = 20, height = 15, units = "cm", dpi = 300)



test_df <- melt(GC06_hdist_present, c("lon", "lat"), value.name = "H_dist")

p6 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon-180, y=lat-90, fill=H_dist))+
  labs(subtitle = 'Smooth = 0.6, Calibration period : 1950 - 1975')+
  ggtitle(paste0('GC BC', ': Average H = ', round(mean(GC06_hdist_present), 2)))+
  scale_fill_gradient(low = "white", high = "#015a8c", limits = c(0.1, 0.70), oob = scales::squish)+
  borders("world", colour = 'black', lwd = 0.12) +
  scale_x_continuous(, expand = c(0, 0)) +
  scale_y_continuous(, expand = c(0, 0))+
  theme(legend.position = 'bottom')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())+
  xlab('Longitude')+
  ylab('Latitude') +
  labs(fill='Hellinger \nDistance')+
  theme_bw()+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=16), #change legend title font sizen
        legend.text = element_text(size=12))+ #change legend text font size
  theme(plot.title = element_text(size=24),
        plot.subtitle = element_text(size = 20,hjust=0.5),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),)+
  easy_center_title()
p6

# Generate file name based on the smooth cost
name <- paste0("figure/GC_H_dist_smooth06_BC_calibration")

# Save the plot as both PDF and PNG
ggsave(paste0(name, ".pdf"), plot = p6, width = 20, height = 15, units = "cm", dpi = 300)
ggsave(paste0(name, ".png"), plot = p6, width = 20, height = 15, units = "cm", dpi = 300)



# Map MMM H dist
#Present
test_df <- melt(MMM_hdist, varnames = c("lon", "lat"), value.name = "H_dist")

p6 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon-180, y=lat-90, fill=H_dist))+
  labs(subtitle = 'Projection period : 1998 - 2023')+
  ggtitle(paste0('MMM', ': Average H = ', round(mean(MMM_hdist), 2)))+
  scale_fill_gradient(low = "white", high = "#015a8c", limits = c(0.1, 0.70), oob = scales::squish)+
  borders("world", colour = 'black', lwd = 0.12) +
  scale_x_continuous(, expand = c(0, 0)) +
  scale_y_continuous(, expand = c(0,0))+
  theme(legend.position = 'bottom')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())+
  xlab('Longitude')+
  ylab('Latitude') +
  labs(fill='Hellinger \nDistance')+
  theme_bw()+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=16), #change legend title font sizen
        legend.text = element_text(size=12))+ #change legend text font size
  theme(plot.title = element_text(size=24),
        plot.subtitle = element_text(size = 20,hjust=0.5),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),)+
  easy_center_title()
p6

# Generate file name based on the smooth cost
name <- paste0("figure/MMM_H_dist_BC_calibration")

# Save the plot as both PDF and PNG
ggsave(paste0(name, ".pdf"), plot = p6, width = 20, height = 15, units = "cm", dpi = 300)
ggsave(paste0(name, ".png"), plot = p6, width = 20, height = 15, units = "cm", dpi = 300)

#Future
# Melt it into a data frame:
test_df <- melt(MMM_hdist_future, varnames = c("lon", "lat"), value.name = "H_dist")

p6 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon-180, y=lat-90, fill=H_dist))+
  labs(subtitle = 'Projection period : 1998 - 2023')+
  ggtitle(paste0('MMM', ': Average H = ', round(mean(MMM_hdist_future), 2)))+
  scale_fill_gradient(low = "white", high = "#015a8c", limits = c(0.1, 0.70), oob = scales::squish)+
  borders("world", colour = 'black', lwd = 0.12) +
  scale_x_continuous(, expand = c(0, 0)) +
  scale_y_continuous(, expand = c(0,0))+
  theme(legend.position = 'bottom')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())+
  xlab('Longitude')+
  ylab('Latitude') +
  labs(fill='Hellinger \nDistance')+
  theme_bw()+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=16), #change legend title font sizen
        legend.text = element_text(size=12))+ #change legend text font size
  theme(plot.title = element_text(size=24),
        plot.subtitle = element_text(size = 20,hjust=0.5),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),)+
  easy_center_title()
p6

# Generate file name based on the smooth cost
name <- paste0("figure/MMM_H_dist_BC_projection")

# Save the plot as both PDF and PNG
ggsave(paste0(name, ".pdf"), plot = p6, width = 20, height = 15, units = "cm", dpi = 300)
ggsave(paste0(name, ".png"), plot = p6, width = 20, height = 15, units = "cm", dpi = 300)




# Map H Dist future


test_df <- melt(h_dist_future[,,1], c("lon", "lat"), value.name = "H_dist")

p6 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon-180, y=lat-90, fill=H_dist))+
  labs(subtitle = 'Projection period : 1998 - 2023')+
  ggtitle(paste0('GC Bias corrected', ': Average Hellinger distance = ', round(mean(h_dist_future[,,1]), 2)))+
  scale_fill_gradient(low = "white", high = "#015a8c", limits = c(0.1, 0.70), oob = scales::squish)+
  scale_x_continuous(, expand = c(0, 0)) +
  scale_y_continuous(, expand = c(0, 0))+
  theme(legend.position = 'bottom')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())+
  xlab('Longitude')+
  ylab('Latitude') +
  labs(fill='Hellinger \nDistance')+
  theme_bw()+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=16), #change legend title font sizen
        legend.text = element_text(size=12))+ #change legend text font size
  theme(plot.title = element_text(size=24),
        plot.subtitle = element_text(size = 20,hjust=0.5),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),)+
  easy_center_title()
p6



{

  for (m in seq_along(model_names)) {
    model_hdist <- h_dist_future[, , m]

    # Melt the matrix into a dataframe with lon/lat
    test_df <- melt(model_hdist, varnames = c("lon_idx", "lat_idx"), value.name = "H_dist")
    test_df$lon <- lon[test_df$lon_idx]
    test_df$lat <- lat[test_df$lat_idx]

    # Plot
    p <- ggplot() +
      geom_tile(data = test_df, aes(x = lon, y = lat, fill = H_dist)) +
      labs(subtitle = 'Projection period : 1998 - 2023') +
      ggtitle(paste0(model_names[m], ': Average H = ',
                     round(mean(model_hdist, na.rm = TRUE), 2))) +
      scale_fill_gradient(low = "white", high = "#015a8c", limits = c(0.1, 0.70), oob = scales::squish) +
      borders("world", colour = 'black', lwd = 0.12) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      theme(legend.position = 'bottom') +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank()) +
      xlab('Longitude') +
      ylab('Latitude') +
      labs(fill = 'Hellinger \nDistance') +
      theme_bw() +
      theme(legend.key.size = unit(1, 'cm'),
            legend.key.height = unit(1.4, 'cm'),
            legend.key.width = unit(0.4, 'cm'),
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 12),
            plot.title = element_text(size = 24),
            plot.subtitle = element_text(size = 20, hjust = 0.5),
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 16)) +
      easy_center_title()

    print(p)
  }
}

# Map of H dist for each model
{

  for (m in seq_along(model_names)) {
    model_hdist <- partial_hdist_future[, , m]

    # Melt the matrix into a dataframe with lon/lat
    test_df <- melt(model_hdist, varnames = c("lon_idx", "lat_idx"), value.name = "H_dist")
    test_df$lon <- lon[test_df$lon_idx]
    test_df$lat <- lat[test_df$lat_idx]

    # Plot
    p <- ggplot() +
      geom_tile(data = test_df, aes(x = lon, y = lat, fill = H_dist)) +
      labs(subtitle = 'Projection period : 1998 - 2023') +
      ggtitle(paste0(model_names[m], ': Average partial (0.10) H = ',
                     round(mean(model_hdist, na.rm = TRUE), 2))) +
      scale_fill_gradient(low = "white", high = "#015a8c", limits = c(0.1, 0.3), oob = scales::squish) +
      borders("world", colour = 'black', lwd = 0.12) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      theme(legend.position = 'bottom') +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank()) +
      xlab('Longitude') +
      ylab('Latitude') +
      labs(fill = 'Hellinger \nDistance') +
      theme_bw() +
      theme(legend.key.size = unit(1, 'cm'),
            legend.key.height = unit(1.4, 'cm'),
            legend.key.width = unit(0.4, 'cm'),
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 12),
            plot.title = element_text(size = 24),
            plot.subtitle = element_text(size = 20, hjust = 0.5),
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 16)) +
      easy_center_title()

    print(p)
  }
}

# Inspect psl in era5
{
  library(ncdf4)

  # User input
  lon <- 0:180
  lat <- -90:90
  lon_idx <- 9
  lat_idx <- 7
  year_range <- 1950:1975
  var <- "psl"

  # Open NetCDF
  f <- "data/CMIP6_merged_all/CMCC-ESM2/psl/psl_CMCC-ESM2_19500101-21001230.nc"
  nc <- nc_open(f)

  # Get grid
  lon_file <- ncvar_get(nc, "lon")
  lat_file <- ncvar_get(nc, "lat")
  lon_file_adjusted <- ifelse(lon_file > 180, lon_file - 360, lon_file)
  lon_user_adjusted <- ifelse(lon > 180, lon - 360, lon)

  # Match user-defined region to file grid indices
  lon_order <- order(lon_file_adjusted)
  lon_file_sorted <- lon_file_adjusted[lon_order]
  lon_file_original_sorted <- lon_file[lon_order]
  lon_idx_unsorted <- match(lon_user_adjusted, lon_file_sorted)
  lon_idx_in_file <- lon_order[lon_idx_unsorted]
  lat_idx_in_file <- match(lat, lat_file)

  # Extract date and match time indices
  yyyy <- extract_years_from_time(nc)  # assuming this function is already defined
  iyears <- which(yyyy %in% year_range)

  # Compute final indices
  start_lon <- min(lon_idx_in_file)
  count_lon <- max(lon_idx_in_file) - start_lon + 1
  start_lat <- min(lat_idx_in_file)
  count_lat <- max(lat_idx_in_file) - start_lat + 1
  start_time <- min(iyears)
  count_time <- max(iyears) - start_time + 1

  # Load data block
  data_block <- ncvar_get(nc, var,
                          start = c(start_lon, start_lat, start_time),
                          count = c(count_lon, count_lat, count_time))

  # Map user lon/lat idx to local offset
  lon_file_subset <- lon_file[start_lon:(start_lon + count_lon - 1)]
  lat_file_subset <- lat_file[start_lat:(start_lat + count_lat - 1)]
  lon_idx_local <- match(lon_file_original_sorted[lon_idx_unsorted], lon_file_subset)
  lat_idx_local <- match(lat, lat_file_subset)

  # Extract time series
  ts <- data_block[lon_idx_local[lon_idx], lat_idx_local[lat_idx], ]


  nc_close(nc)

}

# Summary stats
{
  # Dimensions
  nlon_local <- dim(data_block)[1]
  nlat_local <- dim(data_block)[2]

  # Initialize result matrices
  mean_mat <- matrix(NA, nrow = nlon_local, ncol = nlat_local)
  sd_mat   <- matrix(NA, nrow = nlon_local, ncol = nlat_local)
  min_mat  <- matrix(NA, nrow = nlon_local, ncol = nlat_local)
  max_mat  <- matrix(NA, nrow = nlon_local, ncol = nlat_local)

  # Loop over grid points
  for (i in seq_len(nlon_local)) {
    for (j in seq_len(nlat_local)) {
      ts <- ts[!is.na(ts)]  # Remove NAs just in case
      if (length(ts) > 0) {
        mean_mat[i, j] <- mean(data_block[i,j,])
        sd_mat[i, j]   <- sd(data_block[i,j,])
        min_mat[i, j]  <- min(data_block[i,j,])
        max_mat[i, j]  <- max(data_block[i,j,])
      }
    }
  }

}


# Hellinger distance on MV extremes :
{
  # --- Step 1: Compute ldr_indices for each grid point ---
  # The idea is to compute, for each grid point, the indices of the bins that are NOT in the "central" region.
  # The central region is defined using select_hdr_indices(), which returns the bins containing the central (1-tau) mass.
  # Here, tau=0.1 means we keep the central 90% and treat the remaining 10% as "low density" (i.e., extreme) bins.
  #
  # Assuming:
  #   pdf_ref_future has dimensions [lon_size, lat_size, nbins]
  #   nbins is defined, and lon_size = length(lon), lat_size = length(lat)
  ldr_indices <- vector("list", lon_size)
  for (i in seq_len(lon_size)) {
    ldr_indices[[i]] <- vector("list", lat_size)
    for (j in seq_len(lat_size)) {
      # Compute the central region indices using tau=0.1
      central_indices <- select_hdr_indices(pdf_ref_future[i, j, ], tau = 0.10)
      # Then, define the low density indices as those not in the central region
      ldr_indices[[i]][[j]] <- setdiff(seq_len(nbins), central_indices)
    }
  }

  # --- Step 2: Compute partial Hellinger distances on the ldr regions ---
  # We use compute_partial_hdist() to compute, for each grid point and for each model,
  # the Hellinger distance between the reference and model PDFs, but only over the bins defined by the ldr_indices.
  # This returns an array with dimensions [lon_size, lat_size, n_models]
  partial_hdist_future <- compute_partial_hdist(pdf_ref_future, pdf_models_future, ldr_indices)

  # --- Step 3: (Optional) Compute Mean Partial Hellinger Distances per model ---
  # This gives a summary (scalar) for each model.
  n_models <- length(model_names)
  mean_partial_hdist <- sapply(seq_len(n_models), function(m) {
    mean(partial_hdist_future[,, m], na.rm = TRUE)
  })

  # --- Step 4: (Optional) Compute Full Hellinger Distances for Comparison ---
  # For example, if you have previously computed the full Hellinger distances in an array h_dist_future
  mean_h_dist_future <- sapply(seq_len(n_models), function(m) {
    mean(h_dist_future[,, m], na.rm = TRUE)
  })

  # --- Final Outputs ---
  # partial_hdist: array [lon_size, lat_size, n_models] with the partial Hellinger distance (on ldr bins) at each grid point.
  # mean_partial_hdist: vector of length n_models with average partial distances per model.
  # mean_full_h_dist: vector of length n_models with average full Hellinger distances per model.

  # Print results for inspection:
  print("Mean Partial Hellinger Distance per model:")
  print(mean_partial_hdist)
  print("Mean Full Hellinger Distance per model:")
  print(mean_h_dist_future)


  # Ensure model_names is a character vector
  model_names <- as.character(model_names)

  # Create a dataframe for plotting
  df_hdist <- data.frame(
    Model = rep(model_names, 2),
    Mean_Hellinger_Distance = c(mean_h_dist_future, mean_partial_hdist),
    Type = rep(c("Full Hellinger Distance", "Partial Hellinger (Low Density Regions (10%))"), each = length(model_names)),
    stringsAsFactors = FALSE
  )

  # Plot the results
  ggplot(df_hdist, aes(x = Model, y = Mean_Hellinger_Distance, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    theme_minimal() +
    scale_fill_manual(values = c("steelblue", "darkorange")) +
    labs(
      title = "Comparison of Full vs. Partial Hellinger Distance per Model",
      x = "Climate Model",
      y = "Mean Hellinger Distance",
      fill = "Distance Type"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.title = element_text(size = 12),
          legend.position = "bottom")


  # Compute partial H dist on the GC results :
  # H Dist GC
  # Initialize a lon x lat matrix for each smooth cost
  GC01_partial_hdist_present <- matrix(NA, nrow = length(lon), ncol = length(lat))
  GC01_partial_hdist_future <- matrix(NA, nrow = length(lon), ncol = length(lat))

  GC06_partial_hdist_present <- matrix(NA, nrow = length(lon), ncol = length(lat))
  GC06_partial_hdist_future <- matrix(NA, nrow = length(lon), ncol = length(lat))

  for(l in 1:(length(model_names))){  # Ensure that indexing aligns with model names
    islabel <- which(GC_result01$label_attribution == l)
    # GC01_partial_hdist_present[islabel] <- h_dist_present[,,l][islabel]
    GC01_partial_hdist_future[islabel] <- partial_hdist_future[,,l][islabel]

    islabel <- which(GC_result06$label_attribution == l)
    # GC06_partial_hdist_present[islabel] <- h_dist_present[,,l][islabel]
    GC06_partial_hdist_future[islabel] <- partial_hdist_future[,,l][islabel]
  }


  # --- Compute Partial Hellinger Distance for MMM (Future) on Low Density Regions ---

  # Initialize the output matrix for partial Hellinger distances (MMM)

  MMM_partial_hdist <- array(NA, dim = c(lon_size, lat_size))

  for (i in seq_len(lon_size)) {
    for (j in seq_len(lat_size)) {
      # Extract the reference PDF vector at grid point (i, j)
      pdf_ref_vec <- pdf_ref_future[i, j, ]
      # Get the selected low density indices for this grid point
      selected_bins <- ldr_indices[[i]][[j]]

      if (length(selected_bins) > 0) {
        MMM_partial_hdist[i, j] <- sqrt(
          sum(( sqrt(MMM_future[i, j, selected_bins]) - sqrt(pdf_ref_vec[selected_bins]) )^2)
        ) / sqrt(2)
      } else {
        MMM_partial_hdist[i, j] <- NA
      }
    }
  }
  # Print summary statistics
  cat("Mean Partial Hellinger Distance for MMM (future):", mean(MMM_partial_hdist, na.rm = TRUE), "\n")





  # MMM Map of partial hellinger distance future
  # Melt it into a data frame:
  test_df <- melt(MMM_partial_hdist, varnames = c("lon", "lat"), value.name = "partial_H_dist")

  p6 <- ggplot() +
    geom_tile(data=test_df, aes(x=lon-180, y=lat-90, fill=partial_H_dist))+
    labs(subtitle = 'Projection period : 1998 - 2023')+
    ggtitle(paste0('MMM', ': Avg partial H (0.10 LDR) = ', round(mean(MMM_partial_hdist), 2)))+
    scale_fill_gradient(low = "white", high = "#015a8c", limits = c(0.1, 0.3), oob = scales::squish)+
    borders("world", colour = 'black', lwd = 0.12) +
    scale_x_continuous(, expand = c(0, 0)) +
    scale_y_continuous(, expand = c(0,0))+
    theme(legend.position = 'bottom')+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(panel.background = element_blank())+
    xlab('Longitude')+
    ylab('Latitude') +
    labs(fill='Hellinger \nDistance')+
    theme_bw()+
    theme(legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1.4, 'cm'), #change legend key height
          legend.key.width = unit(0.4, 'cm'), #change legend key width
          legend.title = element_text(size=16), #change legend title font sizen
          legend.text = element_text(size=12))+ #change legend text font size
    theme(plot.title = element_text(size=24),
          plot.subtitle = element_text(size = 20,hjust=0.5),
          axis.text=element_text(size=14),
          axis.title=element_text(size=16),)+
    easy_center_title()
  p6

  # Generate file name based on the smooth cost
  name <- paste0("figure/MMM_partial_H_dist_BC")

  # Save the plot as both PDF and PNG
  ggsave(paste0(name, ".pdf"), plot = p6, width = 20, height = 15, units = "cm", dpi = 300)
  ggsave(paste0(name, ".png"), plot = p6, width = 20, height = 15, units = "cm", dpi = 300)


  # GC 01 Map of partial hellinger distance future


  # Melt it into a data frame:
  test_df <- melt(GC01_partial_hdist_future, varnames = c("lon", "lat"), value.name = "partial_H_dist")

  p6 <- ggplot() +
    geom_tile(data=test_df, aes(x=lon-180, y=lat-90, fill=partial_H_dist))+
    labs(subtitle = 'Projection period : 1998 - 2023')+
    ggtitle(paste0('GC lambda = 0.1', ': Avg partial H (0.10 LDR) = ', round(mean(GC01_partial_hdist_future), 2)))+
    scale_fill_gradient(low = "white", high = "#015a8c", limits = c(0.1, 0.3), oob = scales::squish)+
    borders("world", colour = 'black', lwd = 0.12) +
    scale_x_continuous(, expand = c(0, 0)) +
    scale_y_continuous(, expand = c(0,0))+
    theme(legend.position = 'bottom')+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(panel.background = element_blank())+
    xlab('Longitude')+
    ylab('Latitude') +
    labs(fill='Hellinger \nDistance')+
    theme_bw()+
    theme(legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1.4, 'cm'), #change legend key height
          legend.key.width = unit(0.4, 'cm'), #change legend key width
          legend.title = element_text(size=16), #change legend title font sizen
          legend.text = element_text(size=12))+ #change legend text font size
    theme(plot.title = element_text(size=24),
          plot.subtitle = element_text(size = 20,hjust=0.5),
          axis.text=element_text(size=14),
          axis.title=element_text(size=16),)+
    easy_center_title()
  p6

  # Generate file name based on the smooth cost
  name <- paste0("figure/GC01_partial_H_dist_BC")

  # Save the plot as both PDF and PNG
  ggsave(paste0(name, ".pdf"), plot = p6, width = 20, height = 15, units = "cm", dpi = 300)
  ggsave(paste0(name, ".png"), plot = p6, width = 20, height = 15, units = "cm", dpi = 300)


  # GC 06 Map of partial hellinger distance future


  # Melt it into a data frame:
  test_df <- melt(GC06_partial_hdist_future, varnames = c("lon", "lat"), value.name = "partial_H_dist")

  p6 <- ggplot() +
    geom_tile(data=test_df, aes(x=lon-180, y=lat-90, fill=partial_H_dist))+
    labs(subtitle = 'Projection period : 1998 - 2023')+
    ggtitle(paste0('GC lambda = 0.6', ': Avg partial H (0.10 LDR) = ', round(mean(GC06_partial_hdist_future), 2)))+
    scale_fill_gradient(low = "white", high = "#015a8c", limits = c(0.1, 0.3), oob = scales::squish)+
    borders("world", colour = 'black', lwd = 0.12) +
    scale_x_continuous(, expand = c(0, 0)) +
    scale_y_continuous(, expand = c(0,0))+
    theme(legend.position = 'bottom')+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(panel.background = element_blank())+
    xlab('Longitude')+
    ylab('Latitude') +
    labs(fill='Hellinger \nDistance')+
    theme_bw()+
    theme(legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1.4, 'cm'), #change legend key height
          legend.key.width = unit(0.4, 'cm'), #change legend key width
          legend.title = element_text(size=16), #change legend title font sizen
          legend.text = element_text(size=12))+ #change legend text font size
    theme(plot.title = element_text(size=24),
          plot.subtitle = element_text(size = 20,hjust=0.5),
          axis.text=element_text(size=14),
          axis.title=element_text(size=16),)+
    easy_center_title()
  p6

  # Generate file name based on the smooth cost
  name <- paste0("figure/GC06_partial_H_dist_BC")

  # Save the plot as both PDF and PNG
  ggsave(paste0(name, ".pdf"), plot = p6, width = 20, height = 15, units = "cm", dpi = 300)
  ggsave(paste0(name, ".png"), plot = p6, width = 20, height = 15, units = "cm", dpi = 300)

}


# Crossplot of MMM vs GC H dist
{
  # Create a data frame from the grid
  df_cross <- expand.grid(lon = lon, lat = lat)

  # Add the MMM and GC Hellinger distances (flatten the matrices)
  df_cross$MMM <- as.vector(MMM_hdist_future)
  df_cross$GC  <- as.vector(GC01_hdist_future)

  # Clean data (remove rows with NA)
  df_clean <- df_cross[complete.cases(df_cross$MMM, df_cross$GC), ]

  # Compute axis limits from the data
  data_min <- min(df_clean$MMM, df_clean$GC, na.rm = TRUE)
  data_max <- max(df_clean$MMM, df_clean$GC, na.rm = TRUE)
  buffer   <- 0.05 * (data_max - data_min)  # 5% of the range

  min_val <- max(0, data_min - buffer)  # prevent negative if working with distances
  max_val <- data_max

  # Load ggplot2
  library(ggplot2)

  # Plot the crossplot; color the points according to the latitude
  p <- ggplot(df_clean, aes(x = MMM, y = GC, color = lat)) +
    geom_point(alpha = 0.5, size = 0.3) +
    annotate(
      "segment",
      x = min_val, y = min_val, xend = max_val, yend = max_val,
      linetype = "dashed", color = "black", linewidth = 0.5
    ) +
    scale_x_continuous(limits = c(min_val, max_val), expand = c(0, 0)) +
    scale_y_continuous(limits = c(min_val, max_val), expand = c(0, 0)) +
    scale_color_gradient(low = "blue", high = "red") +
    labs(
      title = "Hellinger Distance : MMM | GC lambda = 0.1",
      x = "MMM H",
      y = "GraphCut H",
      color = "Latitude"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 14),
      axis.text  = element_text(size = 12),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )

  p
}

{

  # Create a data frame from the grid
  df_cross <- expand.grid(lon = lon, lat = lat)

  # Add the MMM and GC Hellinger distances (flatten the matrices)
  df_cross$MMM <- as.vector(MMM_partial_hdist)
  df_cross$GC  <- as.vector(GC01_partial_hdist_future)

  # Clean data (remove rows with NA)
  df_clean <- df_cross[complete.cases(df_cross$MMM, df_cross$GC), ]

  # Compute axis limits from the data
  data_min <- min(df_clean$MMM, df_clean$GC, na.rm = TRUE)
  data_max <- max(df_clean$MMM, df_clean$GC, na.rm = TRUE)
  buffer   <- 0.05 * (data_max - data_min)  # 5% of the range

  min_val <- max(0, data_min - buffer)  # prevent negative if working with distances
  max_val <- data_max

  # Load ggplot2
  library(ggplot2)

  # Plot the crossplot; color the points according to the latitude
  ggplot(df_clean, aes(x = MMM, y = GC, color = lat)) +
    geom_point(alpha = 0.5, size = 0.3) +
    annotate(
      "segment",
      x = min_val, y = min_val, xend = max_val, yend = max_val,
      linetype = "dashed", color = "black", linewidth = 0.5
    ) +
    scale_x_continuous(limits = c(min_val, max_val), expand = c(0, 0)) +
    scale_y_continuous(limits = c(min_val, max_val), expand = c(0, 0)) +
    scale_color_gradient(low = "blue", high = "red") +
    labs(
      title = "Partial H (LDR Mass = 0.1) : MMM | GC lambda = 0.1",
      x = "MMM H",
      y = "GraphCut H",
      color = "Latitude"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 14),
      axis.text  = element_text(size = 12),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )

}


# Crossplot of the H dist vs lat
{

  df_cross$MMM <- as.vector(MMM_hdist_future)
  df_cross$GC  <- as.vector(GC01_hdist_future)
  # Clean data
  df_clean <- df_cross[complete.cases(df_cross$MMM, df_cross$GC), ]

  # Convert to long format for plotting
  library(tidyr)
  df_long <- pivot_longer(
    df_clean,
    cols = c(MMM, GC),
    names_to = "Method",
    values_to = "Hellinger"
  )

  # Plot
  library(ggplot2)

  ggplot(df_long, aes(x = Hellinger, y = lat, color = Method)) +
    geom_point(alpha = 0.4, size = 0.3) +
    facet_wrap(~Method, nrow = 1) +
    scale_color_manual(values = c("MMM" = "steelblue", "GC" = "tomato")) +
    labs(
      title = "Hellinger Distance vs Latitude",
      y = "Latitude",
      x = "Hellinger Distance"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 14),
      axis.text  = element_text(size = 12),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      strip.text = element_text(size = 13)
    )

}


# Crossplot of the diff in H between MMM and GC by latitude with median
{
  diff_h <- MMM_hdist_future - GC01_hdist_future

  # Create a data frame of all lon/lat pairs
  df_test <- expand.grid(lon = lon, lat = lat)

  # Flatten 'test' and assign to dataframe
  df_test$Hellinger <- as.vector(diff_h)

  # Remove missing values
  df_test <- df_test[complete.cases(df_test$Hellinger), ]

  # Add label for better method
  df_test$Better <- ifelse(df_test$Hellinger < 0, "MMM better", "GC better")

  # Compute median diff(H) per latitude
  library(dplyr)
  median_lat <- df_test %>%
    group_by(lat) %>%
    summarise(Hellinger = median(Hellinger, na.rm = TRUE)) %>%
    mutate(Better = "Median")

  # Combine data for consistent color mapping
  df_plot <- bind_rows(df_test, median_lat)

  # Plot
  library(ggplot2)

  ggplot(df_plot, aes(x = Hellinger, y = lat, color = Better)) +
    geom_point(data = df_test, alpha = 1, size = 0.3) +
    geom_point(data = median_lat, size = 0.8) +
    scale_color_manual(
      values = c("MMM better" = "blue", "GC better" = "red", "Median" = "black")
    ) +
    labs(
      title = "Difference of H by latitude (H(MMM) - H(GC))",
      y = "Latitude",
      x = "diff(H)"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 14),
      axis.text  = element_text(size = 12),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 13),
      legend.text  = element_text(size = 12)
    )
}

# Crossplot of the diff in H between MMM and GC by latitude with median
{
  diff_h <- MMM_partial_hdist - GC01_partial_hdist_future

  # Create a data frame of all lon/lat pairs
  df_test <- expand.grid(lon = lon, lat = lat)

  # Flatten 'test' and assign to dataframe
  df_test$Hellinger <- as.vector(diff_h)

  # Remove missing values
  df_test <- df_test[complete.cases(df_test$Hellinger), ]

  # Add label for better method
  df_test$Better <- ifelse(df_test$Hellinger < 0, "MMM better", "GC better")

  # Compute median diff(H) per latitude
  library(dplyr)
  median_lat <- df_test %>%
    group_by(lat) %>%
    summarise(Hellinger = median(Hellinger, na.rm = TRUE)) %>%
    mutate(Better = "Median")

  # Combine data for consistent color mapping
  df_plot <- bind_rows(df_test, median_lat)

  # Plot
  library(ggplot2)

  ggplot(df_plot, aes(x = Hellinger, y = lat, color = Better)) +
    geom_point(data = df_test, alpha = 1, size = 0.3) +
    geom_point(data = median_lat, size = 0.8) +
    scale_color_manual(
      values = c("MMM better" = "blue", "GC better" = "red", "Median" = "black")
    ) +
    labs(
      title = "Difference of partial H by latitude (H(MMM) - H(GC))",
      y = "Latitude",
      x = "diff(H)"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 14),
      axis.text  = element_text(size = 12),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 13),
      legend.text  = element_text(size = 12)
    )
}


# Crossplot of the diff in H between MMM and GC by longitude with median
{
  diff_h <- MMM_hdist_future - GC01_hdist_future
  # Create a data frame of all lon/lat pairs
  df_test <- expand.grid(lon = lon, lat = lat)

  # Flatten 'test' and assign to dataframe
  df_test$Hellinger <- as.vector(diff_h)

  # Remove missing values
  df_test <- df_test[complete.cases(df_test$Hellinger), ]

  # Add label for better method
  df_test$Better <- ifelse(df_test$Hellinger < 0, "MMM better", "GC better")

  # Compute median diff(H) per longitude
  library(dplyr)
  median_lon <- df_test %>%
    group_by(lon) %>%
    summarise(Hellinger = median(Hellinger, na.rm = TRUE)) %>%
    mutate(Better = "Median")

  # Combine data for consistent color mapping
  df_plot <- bind_rows(df_test, median_lon)

  # Plot
  library(ggplot2)

  ggplot(df_plot, aes(x = lon, y = Hellinger, color = Better)) +
    geom_point(data = df_test, alpha = 1, size = 0.3) +
    geom_point(data = median_lon, size = 0.8) +
    scale_color_manual(
      values = c("MMM better" = "blue", "GC better" = "red", "Median" = "black")
    ) +
    labs(
      title = "Difference of H by longitude (H(MMM) - H(GC))",
      x = "Longitude",
      y = "diff(H)"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 14),
      axis.text  = element_text(size = 12),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 13),
      legend.text  = element_text(size = 12)
    )
}


# Crossplot of the diff in H between MMM and GC by longitude with median
{
  diff_h <- MMM_partial_hdist - GC01_partial_hdist_future
  # Create a data frame of all lon/lat pairs
  df_test <- expand.grid(lon = lon, lat = lat)

  # Flatten 'test' and assign to dataframe
  df_test$Hellinger <- as.vector(diff_h)

  # Remove missing values
  df_test <- df_test[complete.cases(df_test$Hellinger), ]

  # Add label for better method
  df_test$Better <- ifelse(df_test$Hellinger < 0, "MMM better", "GC better")

  # Compute median diff(H) per longitude
  library(dplyr)
  median_lon <- df_test %>%
    group_by(lon) %>%
    summarise(Hellinger = median(Hellinger, na.rm = TRUE)) %>%
    mutate(Better = "Median")

  # Combine data for consistent color mapping
  df_plot <- bind_rows(df_test, median_lon)

  # Plot
  library(ggplot2)

  ggplot(df_plot, aes(x = lon, y = Hellinger, color = Better)) +
    geom_point(data = df_test, alpha = 1, size = 0.3) +
    geom_point(data = median_lon, size = 0.8) +
    scale_color_manual(
      values = c("MMM better" = "blue", "GC better" = "red", "Median" = "black")
    ) +
    labs(
      title = "Difference of partial H by longitude (H(MMM) - H(GC))",
      x = "Longitude",
      y = "diff(H)"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 14),
      axis.text  = element_text(size = 12),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 13),
      legend.text  = element_text(size = 12)
    )
}