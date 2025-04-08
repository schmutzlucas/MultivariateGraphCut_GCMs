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

test_hdist <- compute_partial_hdist(pdf_ref_future[180 + 131, 90 -7, ], pdf_ref_present[180 + 131, 90 -7, ], 1:512)

# Plotting the 3D PDF of one gridpoint
{
  library(plotly)
  library(abind)

  # Example indices (adjust as needed)
  lon_index <- 2  # selected longitude index
  lat_index <- 2  # selected latitude index
  model_idx <- 2  # select one model (from the pdf_models_future array)
  nbins <- 8      # number of bins per variable

  # Extract the PDF vector for the chosen grid point and model from the results.
  pdf_vector <- results$pdf_ref_present[lon_index, lat_index, ]

  # Reshape the 1D PDF vector into a 3D array.
  pdf_3d <- array(pdf_vector, dim = c(nbins, nbins, nbins))

  # Adjust range_var to extract the correct reference range for the chosen grid point
  range_var_adj <- aperm(abind(range_var_final$ranges, along = 4), c(1, 2, 4, 3))
  range_mat <- range_var_adj[lon_index + 130, lat_index + 90, , ]  # dimensions: [n_vars, 2]

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
  lon_index <- 1 # selected longitude index
  lat_index <- 3 # selected latitude index
  model_idx <- 1 # select one model (from the pdf_models array)
  nbins <- 8       # number of bins per variable

  # Extract the PDF vector for the chosen grid point and model.
  # Here pdf_models is from results$pdf_models$present and has dimensions:
  # [lon, lat, nbins^n_vars, num_models].
  pdf_vector <- results$pdf_ref$future[lon_index, lat_index, ]

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

# Total number of bins (joint PDF) per grid point
n_bins_total <- nbins1d^(length(variables))

# Create a nested list of selected indices for every grid point (using all bins)
all_bins <- 1:n_bins_total
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



{
  # Map of labels
  # Extract the label attribution for the current smooth cost
  GC_labels <- GC_result06$label_attribution

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
    ggtitle(paste("Label GC Hellinger - Lambda:", 0.6)) +
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
name <- paste0("figure/GC_labelling_smooth01_BC_present_22model")

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
GC_hdist_present <- matrix(NA, nrow = length(lon), ncol = length(lat))
GC_hdist_future <- matrix(NA, nrow = length(lon), ncol = length(lat))

for(l in 1:(length(model_names))){  # Ensure that indexing aligns with model names
  islabel <- which(GC_result01$label_attribution == l)
  GC_hdist_present[islabel] <- h_dist_present[,,l][islabel]
  GC_hdist_future[islabel] <- h_dist_future[,,l][islabel]
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
    MMM_hdist[i, j] <- sqrt(sum((sqrt(MMM_present[i, j, ]) - sqrt(pdf_ref_present[i, j, ]))^2)) / sqrt(2)

    # Compute Hellinger distance for future
    MMM_hdist_future[i, j] <- sqrt(sum((sqrt(MMM_future[i, j, ]) - sqrt(pdf_ref_future[i, j, ]))^2)) / sqrt(2)
  }
}

# Get the current date and time
current_time <- Sys.time()

# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")

# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "_my_workspace_ERA5_bias_corrected_7models_final.RData")

# Save the workspace using the generated filename
save.image(file = filename, compress = FALSE)



# Map H Dist future


test_df <- melt(GC_hdist_future, c("lon", "lat"), value.name = "H_dist")

p6 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon-180, y=lat-90, fill=H_dist))+
  labs(subtitle = 'Smooth = 0.6, Projection period : 1998 - 2023')+
  ggtitle(paste0('GC Bias corrected', ': Average H = ', round(mean(GC_hdist_future), 2)))+
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
name <- paste0("figure/GC_H_dist_smooth06_BC")

# Save the plot as both PDF and PNG
ggsave(paste0(name, ".pdf"), plot = p6, width = 20, height = 15, units = "cm", dpi = 300)
ggsave(paste0(name, ".png"), plot = p6, width = 20, height = 15, units = "cm", dpi = 300)


# Map H Dist present


test_df <- melt(GC_hdist_present, c("lon", "lat"), value.name = "H_dist")

p6 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon-180, y=lat-90, fill=H_dist))+
  labs(subtitle = 'Smooth = 0.6, Calibration period : 1950 - 1975')+
  ggtitle(paste0('GC BC', ': Average H = ', round(mean(GC_hdist_present), 2)))+
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
name <- paste0("figure/GC_H_dist_smooth06_BC_present")

# Save the plot as both PDF and PNG
ggsave(paste0(name, ".pdf"), plot = p6, width = 20, height = 15, units = "cm", dpi = 300)
ggsave(paste0(name, ".png"), plot = p6, width = 20, height = 15, units = "cm", dpi = 300)



# Map MMM H dist


# Melt it into a data frame:
test_df <- melt(MMM_hdist_future, varnames = c("lon", "lat"), value.name = "H_dist")

p6 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon-180, y=lat-90, fill=H_dist))+
  labs(subtitle = 'Projection period : 1998 - 2023')+
  ggtitle(paste0('MMM', ': Average Hellinger distance = ', round(mean(MMM_hdist_future), 2)))+
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
name <- paste0("figure/MMM_H_dist_BC")

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
  
}