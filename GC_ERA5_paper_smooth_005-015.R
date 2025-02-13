# Install and load necessary libraries
list_of_packages <- read.table("package_list.txt", sep="\n")$V1
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages))
  install.packages(new.packages, repos = "https://cloud.r-project.org")

library(devtools)
lapply(list_of_packages, library, character.only = TRUE)
install_github("schmutzlucas/gcoWrapR")


# Loading local functions
source_code_dir <- 'functions/' #The directory where all functions are saved.
file_paths <- list.files(source_code_dir, full.names = T)
for(path in file_paths){source(path)}

range_var_final <- readRDS('ranges/range_var_final_GreenwichCentered_1950-2023_90deg_3v.rds')

# Setting global variables
lon <- -180:179
lat <- -90:90
lon_size <- length(lon)
lat_size <- length(lat)
# Temporal ranges
year_present <<- 1950:1975
year_future <<- 1998:2023
# data directory
data_dir <<- 'data/CMIP6_merged_all/'

# List of the variable used
variables <- c('pr', 'tas', 'psl')

# Bins for the pdfs
nbins1d <<- 8
nbins <<- nbins1d^(length(variables))


# Obtains the list of models from the model names or from a file
model_names <- read.table('model_names_pr_tas_psl.txt')
model_names <- as.list(model_names[['V1']])
# Index of the reference
ref_index <<- 1
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

# Time the execution of the optimized function
time_optimized <- system.time({
  # todo add number of workers as argument
  tmp <- compute_nd_pdf_optimized_0centered(variables, model_names, data_dir, year_present, year_future,
                                  lon, lat, aperm(abind(range_var_final$ranges, along = 4), c(1, 2, 4, 3)), nbins1d, workers = 3)
})
cat("Time taken for compute_nd_pdf_optimized: ", format_time(time_optimized["elapsed"]), "\n")

# Choose the reference in the models
reference_name <<- model_names[ref_index]
model_names <<- model_names[-ref_index]


# Get the current date and time
current_time <- Sys.time()

# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")

# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "_my_workspace_ERA5_allModels_beforeOptim_3v.RData")

# Save the workspace using the generated filename
save.image(file = filename, compress = FALSE)

pdf_present <- tmp$present
pdf_future <- tmp$future

pdf_ref_present <- pdf_present[ , , , ref_index]
pdf_models_present <- pdf_present[ , , , -ref_index]

pdf_ref_future <- pdf_future[ , , , ref_index]
pdf_models_future <- pdf_future[ , , , -ref_index]

rm(pdf_present, pdf_future)


# Computing the sum of hellinger distances between models and reference --> used as datacost
h_dist <- array(data = 0, dim = c(length(lon), length(lat),
                                  length(model_names)))
h_dist_unchecked <- array(data = 0, dim = c(length(lon), length(lat),
                                            length(model_names)))

h_dist_future <- array(data = 0, dim = c(length(lon), length(lat),
                                         length(model_names)))
h_dist_unchecked_future <- array(data = 0, dim = c(length(lon), length(lat),
                                                   length(model_names)))

# Loop through variables and models
m <- 1
for (model_name in model_names) {
  for (i in seq_along(lon)) {
    for (j in seq_along(lat)) {
      # Compute Hellinger distance
      h_dist_unchecked[i, j, m] <- sqrt(sum((sqrt(pdf_models_present[i, j, , m]) - sqrt(pdf_ref_present[i, j, ]))^2)) / sqrt(2)
      h_dist_unchecked_future[i, j, m] <- sqrt(sum((sqrt(pdf_models_future[i, j, , m]) - sqrt(pdf_ref_future[i, j, ]))^2)) / sqrt(2)
    }
  }
  m <- m + 1
}

hist(h_dist_unchecked)
# Replace NaN with 0
h_dist[,,] <- replace(h_dist_unchecked[,,], is.nan(h_dist_unchecked), 0)
h_dist_future[,,] <- replace(h_dist_unchecked_future[,,], is.nan(h_dist_unchecked), 0)
hist(h_dist)
rm(h_dist_unchecked, h_dist_unchecked_future)


# Get the current date and time
current_time <- Sys.time()

# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")

# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "_my_workspace_ERA5_allModels_beforeOptim_3v.RData")

# Save the workspace using the generated filename
save.image(file = filename, compress = FALSE)


smooth_cost <- 1
# Wrap each iteration in tryCatch to handle errors gracefully
tryCatch({
  # Run Graph Cut with the varying smooth cost
  GC_result11 <- GraphCutHellinger_nD(
    pdf_models_future = pdf_models_present,
    h_dist = h_dist,
    weight_data = 1,               # Fixed data weight
    weight_smooth = smooth_cost,   # Varying smooth cost
    nBins = nbins1d^3,
    seed = 1,
    verbose = TRUE,
    rebuild = TRUE
  )


}, error = function(e) {
  cat("Error encountered with smooth cost =", smooth_cost, ": ", e$message, "\n")
})


smooth_cost <- 0.05
# Wrap each iteration in tryCatch to handle errors gracefully
tryCatch({
  # Run Graph Cut with the varying smooth cost
  GC_result_0.05 <- GraphCutHellinger_nD(
    pdf_models_future = pdf_models_present,
    h_dist = h_dist,
    weight_data = 1,               # Fixed data weight
    weight_smooth = smooth_cost,   # Varying smooth cost
    nBins = nbins1d^3,
    seed = 1,
    verbose = TRUE,
    rebuild = FALSE
  )


}, error = function(e) {
  cat("Error encountered with smooth cost =", smooth_cost, ": ", e$message, "\n")
})

smooth_cost <- 0.1
# Wrap each iteration in tryCatch to handle errors gracefully
tryCatch({
  # Run Graph Cut with the varying smooth cost
  GC_result_0.1 <- GraphCutHellinger_nD(
    pdf_models_future = pdf_models_present,
    h_dist = h_dist,
    weight_data = 1,               # Fixed data weight
    weight_smooth = smooth_cost,   # Varying smooth cost
    nBins = nbins1d^3,
    seed = 1,
    verbose = TRUE,
    rebuild = FALSE
  )


}, error = function(e) {
  cat("Error encountered with smooth cost =", smooth_cost, ": ", e$message, "\n")
})


smooth_cost <- 0.15
# Wrap each iteration in tryCatch to handle errors gracefully
tryCatch({
  # Run Graph Cut with the varying smooth cost
  GC_result_0.15 <- GraphCutHellinger_nD(
    pdf_models_future = pdf_models_present,
    h_dist = h_dist,
    weight_data = 1,               # Fixed data weight
    weight_smooth = smooth_cost,   # Varying smooth cost
    nBins = nbins1d^3,
    seed = 1,
    verbose = TRUE,
    rebuild = FALSE
  )


}, error = function(e) {
  cat("Error encountered with smooth cost =", smooth_cost, ": ", e$message, "\n")
})

# Initialize a lon x lat matrix for each smooth cost
GC_hdist_0.05 <- matrix(NA, nrow = length(lon), ncol = length(lat))
GC_hdist_future_0.05 <- matrix(NA, nrow = length(lon), ncol = length(lat))

for(l in 1:(length(model_names))){  # Ensure that indexing aligns with model names
  islabel <- which(GC_result_0.05$label_attribution == l)
  GC_hdist_0.05[islabel] <- h_dist[,,l][islabel]
  GC_hdist_future_0.05[islabel] <- h_dist_future[,,l][islabel]
}


# Initialize a lon x lat matrix for each smooth cost
GC_hdist_0.1 <- matrix(NA, nrow = length(lon), ncol = length(lat))
GC_hdist_future_0.1 <- matrix(NA, nrow = length(lon), ncol = length(lat))

for(l in 1:(length(model_names))){  # Ensure that indexing aligns with model names
  islabel <- which(GC_result_0.1$label_attribution == l)
  GC_hdist_0.1[islabel] <- h_dist[,,l][islabel]
  GC_hdist_future_0.1[islabel] <- h_dist_future[,,l][islabel]
}

# Initialize a lon x lat matrix for each smooth cost
GC_hdist_0.15 <- matrix(NA, nrow = length(lon), ncol = length(lat))
GC_hdist_future_0.15 <- matrix(NA, nrow = length(lon), ncol = length(lat))

for(l in 1:(length(model_names))){  # Ensure that indexing aligns with model names
  islabel <- which(GC_result_0.15$label_attribution == l)
  GC_hdist_0.15[islabel] <- h_dist[,,l][islabel]
  GC_hdist_future_0.15[islabel] <- h_dist_future[,,l][islabel]
}

hist(GC_hdist_future_0.05)
hist(GC_hdist_future_0.1)
hist(GC_hdist_future_0.15)

mean(GC_hdist_future_0.05)
mean(GC_hdist_future_0.1)
mean(GC_hdist_future_0.15)

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

# Replace NaN values with 0
MMM_hdist <- replace(MMM_hdist, is.nan(MMM_hdist), 0)
MMM_hdist_future <- replace(MMM_hdist_future, is.nan(MMM_hdist_future), 0)

# Visualize or analyze the results
hist(MMM_hdist)
hist(GC_hdist)
hist(MMM_hdist_future)
hist(GC_hdist_future)

mean(MMM_hdist_future)
mean(GC_hdist_future)


# Get the current date and time
current_time <- Sys.time()

# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")

# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "_my_workspace_ERA5_allModels_final_results.RData")

# Save the workspace using the generated filename
save.image(file = filename, compress = FALSE)

# Step 1: Compute ldr_indices (outlier bins for each grid point)
ldr_indices <- vector("list", lon_size)

for (i in 1:lon_size) {  # Longitude first (memory efficient)
  ldr_indices[[i]] <- vector("list", lat_size)  # Initialize sublists

  for (j in 1:lat_size) {  # Latitude next
    hdr_indices <- select_hdr_indices(pdf_ref_future[i, j, ], tau = 0.1)
    ldr_indices[[i]][[j]] <- setdiff(seq_len(nbins), hdr_indices)  # Compute outlier bins
  }
}

# Step 2: Compute partial Hellinger distance using precomputed indices
h_dist_partial <- compute_partial_hdist(pdf_ref_future, pdf_models_future, ldr_indices)


# Convert model names to a character vector (if needed)
model_names <- as.character(model_names)

# Compute mean Hellinger distance (full)
mean_h_dist_future <- sapply(seq_along(model_names), function(m) {
  mean(h_dist_future[,,m], na.rm = TRUE)  # Avoid NaNs
})

# Compute mean Hellinger distance (partial - lower 10% mass)
mean_h_dist_partial <- sapply(seq_along(model_names), function(m) {
  mean(h_dist_partial[,,m], na.rm = TRUE)  # Avoid NaNs
})

# Create a dataframe
df_hdist <- data.frame(
  Model = rep(model_names, 2),  # Ensure Model is a character vector
  Mean_Hellinger_Distance = c(mean_h_dist_future, mean_h_dist_partial),
  Type = rep(c("Full Hellinger Distance", "Partial Hellinger (Low Density Regions (10%))"), each = length(model_names))
)

# Plot
ggplot(df_hdist, aes(x = Model, y = Mean_Hellinger_Distance, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  theme_minimal() +
  scale_fill_manual(values = c("steelblue", "darkorange")) +  # Custom colors
  labs(
    title = "Comparison of Full vs. Partial Hellinger Distance per Model",
    x = "Climate Model",
    y = "Mean Hellinger Distance",
    fill = "Distance Type"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title = element_text(size = 12),
        legend.position = "bottom")


# Step 2: Compute partial Hellinger distance using precomputed indices
h_dist_partial_MMM <- compute_partial_hdist(pdf_ref_future, MMM_future, ldr_indices)


# Melt it into a data frame:
test_df <- melt(h_dist_partial_MMM, varnames = c("lon", "lat"), value.name = "H_dist")

p6 <- ggplot() +
  geom_tile(data = test_df, aes(x = lon, y = lat, fill = H_dist)) +
  labs(subtitle = 'Projection period : 1999 - 2014') +
  ggtitle(paste0('MMM', ': Average Hellinger distance = ',
                 round(mean(h_dist_partial_MMM), 2))) +
  scale_fill_gradient(low = "white", high = "#015a8c",
                      limits = c(0.1, 0.70), oob = scales::squish) +
  # Use "world" instead of "world2" so that it matches –180..+180
  borders("world", colour = 'black', fill = NA, size = 0.12) +
  # Force x,y limits to match your new longitude and latitude
  scale_x_continuous(limits = c(-180, 180), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-90, 90), expand = c(0, 0)) +
  # Keep or replace with coord_quickmap()
  coord_fixed(ratio = 1) +
  theme(legend.position = 'bottom',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  xlab('Longitude') +
  ylab('Latitude') +
  labs(fill = 'Hellinger \nDistance') +
  theme_bw() +
  theme(legend.key.size = unit(1, 'cm'),    # legend key size
        legend.key.height = unit(1.4, 'cm'),
        legend.key.width = unit(0.4, 'cm'),
        legend.title = element_text(size=16),
        legend.text = element_text(size=12),
        plot.title = element_text(size=24),
        plot.subtitle = element_text(size = 20, hjust=0.5),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16)) +
  easy_center_title()

p6

# Initialize arrays to store the GC-selected PDFs
pdf_GC_present <- array(NA, dim = c(lon_size, lat_size, nbins))
pdf_GC_future  <- array(NA, dim = c(lon_size, lat_size, nbins))

# Assign the selected PDFs at each grid point
for (i in 1:lon_size) {
  for (j in 1:lat_size) {
    label <- GC_result$label_attribution[i, j]  # Get the assigned model (1-based index)

    if (!is.na(label)) {  # Ensure valid labels
      pdf_GC_present[i, j, ] <- pdf_models_present[i, j, , label]  # Extract the right PDF
      pdf_GC_future[i, j, ]  <- pdf_models_future[i, j, , label]   # Extract future PDF
    }
  }
}



# Step 2: Compute partial Hellinger distance using precomputed indices
h_dist_partial_GC <- compute_partial_hdist(pdf_ref_future, pdf_GC_future, ldr_indices)

test_df <- melt(h_dist_partial_GC, c("lon", "lat"), value.name = "H_dist")

p6 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=H_dist))+
  labs(subtitle = 'Projection period : 1999 - 2014')+
  ggtitle(paste0('MMM', ': Average Hellinger distance = ', round(mean(h_dist_partial_GC), 2)))+
  scale_fill_gradient(low = "white", high = "#015a8c", limits = c(0.1, 0.70), oob = scales::squish)+
  borders("world2", colour = 'black', lwd = 0.12) +
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

# Extract the label attribution for the current smooth cost
GC_labels <- GC_result11$label_attribution

# Convert the label matrix to a data frame for plotting
# Make sure 'GC_labels' has dimensions [lon, lat],
# where lon ∈ [-180..179] and lat ∈ [-90..90].
label_df <- reshape2::melt(GC_labels, varnames = c("lon", "lat"), value.name = "label_attribution")

# Remove this if lat is already -90..+90 (no extra shift needed)
# label_df$lat <- label_df$lat - 90

# Convert label_attribution to a factor with all model names
label_df$label_attribution <- factor(
  label_df$label_attribution,
  levels = seq_along(model_names),
  labels = model_names
)

color_palette <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
  "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
  "#393b79", "#5254a3", "#6b6ecf"
)

p <- ggplot() +
  # Plot tiles at each (lon, lat)
  geom_tile(data = label_df, aes(x = lon, y = lat, fill = label_attribution)) +

  # Custom color palette for the labels
  scale_fill_manual(
    values = color_palette,
    na.value = "white",
    guide = guide_legend(title = "Model Names", ncol = 1)
  ) +

  ggtitle(paste("Label GC Hellinger - Lambda:", smooth_cost)) +

  # Use "world" (not "world2"), which expects longitude ∈ [-180..+180]
  borders("world", colour = 'black', size = 0.12) +

  # Force x,y to match the new [-180..+180] and [-90..+90] coordinate system
  # and lock the aspect ratio (1:1 degrees)
  coord_fixed(xlim = c(-180, 180), ylim = c(-90, 90), ratio = 1) +

  # A clean theme
  theme_bw() +

  # Adjust the legend and axis text, etc.
  theme(
    legend.position = 'bottom',
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

p


# Generate file name based on the smooth cost
name <- paste0("figure/Labels_GC_Hellinger_smooth_1950-1975_3v")

# Save the plot as both PDF and PNG
ggsave(paste0(name, ".pdf"), plot = p, width = 20, height = 15, units = "cm", dpi = 300)
ggsave(paste0(name, ".png"), plot = p, width = 20, height = 15, units = "cm", dpi = 300)


{
  library(plotly)

  # Example grid point indices
  lon_index <- 180+6 # Example longitude index
  lat_index <- 90 + 46  # Example latitude index

  # Extract the 512-bin PDF vector for the specific grid point
  pdf_vector <- pdf_ref_present[lon_index, lat_index, ]

  # Reshape the PDF vector into a 3D array of dimensions [8, 8, 8]
  nbins <- 8
  pdf_3d <- array(pdf_vector, dim = c(nbins, nbins, nbins))

  # Extract the ranges for the variables from range_var_final
  range_var <- range_var_final$ranges
  var1_min <- range_var$pr[lon_index, lat_index, 1]
  var1_max <- range_var$pr[lon_index, lat_index, 2]
  var2_min <- range_var$tas[lon_index, lat_index, 1]
  var2_max <- range_var$tas[lon_index, lat_index, 2]
  var3_min <- range_var$psl[lon_index, lat_index, 1]
  var3_max <- range_var$psl[lon_index, lat_index, 2]

  # Create bin edges for each variable
  x_bins <- seq(var1_min, var1_max, length.out = nbins + 1)
  y_bins <- seq(var2_min, var2_max, length.out = nbins + 1)
  z_bins <- seq(var3_min, var3_max, length.out = nbins + 1)

  # Create the coordinates for the centers of the bins
  x_centers <- (x_bins[-1] + x_bins[-length(x_bins)]) / 2
  y_centers <- (y_bins[-1] + y_bins[-length(y_bins)]) / 2
  z_centers <- (z_bins[-1] + z_bins[-length(z_bins)]) / 2

  # Expand the grid of coordinates
  grid <- expand.grid(x = x_centers, y = y_centers, z = z_centers)

  # Flatten the PDF array into a vector
  pdf_flat <- as.vector(pdf_3d)

  # Combine the coordinates with the PDF values
  plot_data <- data.frame(
    x = grid$x,
    y = grid$y,
    z = grid$z,
    value = pdf_flat
  )

  # Normalize PDF values for marker size
  normalized_pdf <- pdf_flat / max(pdf_flat, na.rm = TRUE)

  # Plot the 3D histogram
  fig <- plot_ly(
    data = plot_data,
    x = ~x,
    y = ~y,
    z = ~z,
    type = "scatter3d",
    mode = "markers",
    marker = list(
      size = ~normalized_pdf * 75,  # Adjust size scaling factor
      color = ~value,
      colorscale = "Viridis",
      showscale = TRUE
    ),
    text = ~paste("PDF Value:", round(value, 4))
  ) %>%
    layout(
      scene = list(
        xaxis = list(title = "Variable 1 (pr)"),
        yaxis = list(title = "Variable 2 (tas)"),
        zaxis = list(title = "Variable 3 (psl)")
      ),
      title = paste("3D PDF for Grid Point (Lon:", lon_index, ", Lat:", lat_index, ")")
    )

  # Show the plot
  fig
}



# Given a PDF for a specific grid point
pdf_gridpoint <- pdf_ref_present[lon_index, lat_index, ]

# Sort the PDF values while keeping track of their original indices
sorted_indices <- order(pdf_gridpoint)
sorted_values <- pdf_gridpoint[sorted_indices]

# Compute the cumulative sum of the sorted values
cumulative_sum <- cumsum(sorted_values)

# Find the indices where the cumulative sum first exceeds 0.1
threshold_index <- which(cumulative_sum >= 0.1)[1]

# Get the indices of the bins contributing to this cumulative sum
selected_indices <- sorted_indices[1:threshold_index]
selected_values <- pdf_gridpoint[selected_indices]

# Output
cat("Indices of bins:", selected_indices, "\n")
cat("Selected bin values:", selected_values, "\n")
cat("Sum of selected values:", sum(selected_values), "\n")



{
  library(plotly)

  # Example grid point indices
  lon_index <- 180 + 6 # Example longitude index
  lat_index <- 90 + 46  # Example latitude index

  # Extract the 512-bin PDF vector for the specific grid point
  pdf_vector <- pdf_ref_present[lon_index, lat_index, ]

  # Reshape the PDF vector into a 3D array of dimensions [8, 8, 8]
  nbins <- 8
  pdf_3d <- array(pdf_vector, dim = c(nbins, nbins, nbins))

  # Extract the ranges for the variables from range_var_final
  range_var <- range_var_final$ranges
  var1_min <- range_var$pr[lon_index, lat_index, 1]
  var1_max <- range_var$pr[lon_index, lat_index, 2]
  var2_min <- range_var$tas[lon_index, lat_index, 1]
  var2_max <- range_var$tas[lon_index, lat_index, 2]
  var3_min <- range_var$psl[lon_index, lat_index, 1]
  var3_max <- range_var$psl[lon_index, lat_index, 2]

  # Create bin edges for each variable
  x_bins <- seq(var1_min, var1_max, length.out = nbins + 1)
  y_bins <- seq(var2_min, var2_max, length.out = nbins + 1)
  z_bins <- seq(var3_min, var3_max, length.out = nbins + 1)

  # Create the coordinates for the centers of the bins
  x_centers <- (x_bins[-1] + x_bins[-length(x_bins)]) / 2
  y_centers <- (y_bins[-1] + y_bins[-length(y_bins)]) / 2
  z_centers <- (z_bins[-1] + z_bins[-length(z_bins)]) / 2

  # Expand the grid of coordinates
  grid <- expand.grid(x = x_centers, y = y_centers, z = z_centers)

  # Flatten the PDF array into a vector
  pdf_flat <- as.vector(pdf_3d)

  # Identify the indices of the bins contributing to the 10% smallest cumulative sum
  sorted_indices <- order(pdf_flat)
  cumulative_sum <- cumsum(pdf_flat[sorted_indices])
  threshold_index <- which(cumulative_sum >= 0.1)[1]
  selected_indices <- sorted_indices[1:threshold_index]

  # Filter data for the selected bins
  highlight_data <- data.frame(
    x = grid$x[selected_indices],
    y = grid$y[selected_indices],
    z = grid$z[selected_indices],
    value = pdf_flat[selected_indices]
  )

  # Normalize PDF values for marker size
  normalized_pdf <- highlight_data$value / max(highlight_data$value, na.rm = TRUE)

  # Plot the 3D histogram with highlighted bins
  fig <- plot_ly(
    data = highlight_data,
    x = ~x,
    y = ~y,
    z = ~z,
    type = "scatter3d",
    mode = "markers",
    marker = list(
      size = ~normalized_pdf * 75,  # Adjust size scaling factor
      color = ~value,
      colorscale = "Viridis",
      showscale = TRUE
    ),
    text = ~paste("PDF Value:", round(value, 4))
  ) %>%
    layout(
      scene = list(
        xaxis = list(title = "Variable 1 (pr)"),
        yaxis = list(title = "Variable 2 (tas)"),
        zaxis = list(title = "Variable 3 (psl)")
      ),
      title = paste("3D PDF (Highlighted 10%) for Grid Point (Lon:", lon_index, ", Lat:", lat_index, ")")
    )

  # Show the plot
  fig
}

test_df <- melt(GC_hdist_future_0.05, c("lon", "lat"), value.name = "H_dist")

p5 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=H_dist))+
  labs(subtitle = 'Projection period : 1999 - 2014')+
  ggtitle(paste0('GraphCut MV', ': Mean Hellinger distance = ', round(mean(GC_hdist_future_0.05), 2)))+
  scale_fill_gradient(low = "white", high = "#015a8c", limits = c(0.1, 0.70), oob = scales::squish)+
  borders("world2", colour = 'black', lwd = 0.12) +
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
p5

# Generate file name based on the smooth cost
name <- paste0("figure/Hdist_GCMV_1950-1975_3v_smooth_005")

# Save the plot as both PDF and PNG
ggsave(paste0(name, ".pdf"), plot = p5, width = 20, height = 15, units = "cm", dpi = 300)
ggsave(paste0(name, ".png"), plot = p5, width = 20, height = 15, units = "cm", dpi = 300)




test_df <- melt(GC_hdist_future_0.1, c("lon", "lat"), value.name = "H_dist")

p5 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=H_dist))+
  labs(subtitle = 'Projection period : 1999 - 2014')+
  ggtitle(paste0('GraphCut MV', ': Mean Hellinger distance = ', round(mean(GC_hdist_future_0.1), 2)))+
  scale_fill_gradient(low = "white", high = "#015a8c", limits = c(0.1, 0.70), oob = scales::squish)+
  borders("world2", colour = 'black', lwd = 0.12) +
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
p5

# Generate file name based on the smooth cost
name <- paste0("figure/Hdist_GCMV_1950-1975_3v_smooth_010")

# Save the plot as both PDF and PNG
ggsave(paste0(name, ".pdf"), plot = p5, width = 20, height = 15, units = "cm", dpi = 300)
ggsave(paste0(name, ".png"), plot = p5, width = 20, height = 15, units = "cm", dpi = 300)




test_df <- melt(GC_hdist_future_0.15, c("lon", "lat"), value.name = "H_dist")

p5 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=H_dist))+
  labs(subtitle = 'Projection period : 1999 - 2014')+
  ggtitle(paste0('GraphCut MV', ': Mean Hellinger distance = ', round(mean(GC_hdist_future_0.15), 2)))+
  scale_fill_gradient(low = "white", high = "#015a8c", limits = c(0.1, 0.70), oob = scales::squish)+
  borders("world2", colour = 'black', lwd = 0.12) +
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
p5

# Generate file name based on the smooth cost
name <- paste0("figure/Hdist_GCMV_1950-1975_3v_smooth_015")

# Save the plot as both PDF and PNG
ggsave(paste0(name, ".pdf"), plot = p5, width = 20, height = 15, units = "cm", dpi = 300)
ggsave(paste0(name, ".png"), plot = p5, width = 20, height = 15, units = "cm", dpi = 300)

test_df <- melt(MMM_hdist_future, c("lon", "lat"), value.name = "H_dist")

p6 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=H_dist))+
  labs(subtitle = 'Projection period : 1999 - 2014')+
  ggtitle(paste0('MMM', ': Average Hellinger distance = ', round(mean(MMM_hdist_future), 2)))+
  scale_fill_gradient(low = "white", high = "#015a8c", limits = c(0.1, 0.70), oob = scales::squish)+
  borders("world2", colour = 'black', lwd = 0.12) +
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
name <- paste0("figure/Hdist_MMM_1950-1975_3v")

# Save the plot as both PDF and PNG
ggsave(paste0(name, ".pdf"), plot = p6, width = 20, height = 15, units = "cm", dpi = 300)
ggsave(paste0(name, ".png"), plot = p6, width = 20, height = 15, units = "cm", dpi = 300)



gradient_MMM <- gradient_hdist(MMM_hdist_future)
gradient_GCMV <- gradient_hdist(GC_hdist_future)


hist(MMM_hdist_future, xlim = c(0, 1), ylim = c(0, 25000), main = "Histogram of MMM Hellinger Distances",
     xlab = "Hellinger Distance", ylab = "Frequency", col = "lightblue", border = "black")
hist(GC_hdist_future, xlim = c(0, 1), ylim = c(0, 25000), main = "Histogram of GC MV Hellinger Distances",
     xlab = "Hellinger Distance", ylab = "Frequency", col = "lightblue", border = "black")




# Marginal :
{
  library(ggplot2)

  # Example grid point indices
  lon_index <- 180+6 # Example longitude index
  lat_index <- 90 + 46  # Example latitude index

  # Extract the 512-bin PDF vector for the specific grid point
  pdf_vector <- pdf_ref_present[lon_index, lat_index, ]

  # Reshape the PDF vector into a 3D array of dimensions [8, 8, 8]
  nbins <- 8
  pdf_3d <- array(pdf_vector, dim = c(nbins, nbins, nbins))

  # Extract the ranges for tas (Variable 2)
  range_var <- range_var_final$ranges
  var2_min <- range_var$tas[lon_index, lat_index, 1]
  var2_max <- range_var$tas[lon_index, lat_index, 2]

  # Create bin edges and centers for tas
  tas_bins <- seq(var2_min, var2_max, length.out = nbins + 1)
  tas_centers <- (tas_bins[-1] + tas_bins[-length(tas_bins)]) / 2  # Midpoints

  # Compute the marginal PDF by summing over pr (Var1) and psl (Var3)
  tas_marginal <- apply(pdf_3d, 2, sum)  # Sum over dimensions 1 and 3

  # Normalize the marginal PDF (optional)
  tas_marginal <- tas_marginal / sum(tas_marginal)

  # Create a data frame for ggplot
  marginal_data <- data.frame(
    tas = tas_centers,
    probability = tas_marginal
  )

  # Plot the marginal distribution
  ggplot(marginal_data, aes(x = tas, y = probability)) +
    geom_bar(stat = "identity", fill = "skyblue", color = "black") +
    labs(
      title = paste("Marginal PDF of tas for Grid Point (Lon:", lon_index, ", Lat:", lat_index, ")"),
      x = "Temperature (tas)",
      y = "Probability Density"
    ) +
    theme_minimal()
}


# Level sets HDR based approach of the outlier selection.
{
  library(plotly)

  # ---------------------------
  # 1. Set up the grid point and data
  # ---------------------------

  # Example grid point indices
  lon_index <- 180 + 6  # Example longitude index
  lat_index <- 90 + 46  # Example latitude index

  # Extract the 512-bin PDF vector for the specific grid point
  pdf_vector <- pdf_ref_future[lon_index, lat_index, ]

  # Reshape the PDF vector into a 3D array of dimensions [8, 8, 8]
  nbins <- 8
  pdf_3d <- array(pdf_vector, dim = c(nbins, nbins, nbins))

  # ---------------------------
  # 2. Set up bin edges and centers for each variable
  # ---------------------------
  # Extract the ranges for the variables from range_var_final
  range_var <- range_var_final$ranges
  var1_min <- range_var$pr[lon_index, lat_index, 1]
  var1_max <- range_var$pr[lon_index, lat_index, 2]
  var2_min <- range_var$tas[lon_index, lat_index, 1]
  var2_max <- range_var$tas[lon_index, lat_index, 2]
  var3_min <- range_var$psl[lon_index, lat_index, 1]
  var3_max <- range_var$psl[lon_index, lat_index, 2]

  # Create bin edges for each variable
  x_bins <- seq(var1_min, var1_max, length.out = nbins + 1)
  y_bins <- seq(var2_min, var2_max, length.out = nbins + 1)
  z_bins <- seq(var3_min, var3_max, length.out = nbins + 1)

  # Compute the centers of the bins
  x_centers <- (x_bins[-1] + x_bins[-length(x_bins)]) / 2
  y_centers <- (y_bins[-1] + y_bins[-length(y_bins)]) / 2
  z_centers <- (z_bins[-1] + z_bins[-length(z_bins)]) / 2

  # Expand the grid of coordinates for plotting
  grid <- expand.grid(x = x_centers, y = y_centers, z = z_centers)

  # Flatten the PDF array into a vector
  pdf_flat <- as.vector(pdf_3d)

  # ---------------------------
  # 3. Define the level set for the central region (HDR)
  # ---------------------------
  # We wish to consider the "central" region that contains 90% of the probability mass,
  # so that the outliers represent the remaining 10%.
  tau <- 0.10         # Fraction of outlier mass
  target_mass <- 1 - tau  # 0.90: the mass in the central region Q(tau)

  # Sort the indices in descending order so that the highest density bins come first
  sorted_indices_desc <- order(pdf_flat, decreasing = TRUE)

  # Compute the cumulative sum of the PDF values for the sorted bins
  cumulative_sum_desc <- cumsum(pdf_flat[sorted_indices_desc])

  # Find the smallest index m such that the cumulative sum is >= target_mass (0.90)
  threshold_index <- which(cumulative_sum_desc >= target_mass)[1]

  # The threshold density kappa is the PDF value at that sorted index
  kappa <- pdf_flat[sorted_indices_desc[threshold_index]]

  # Define the central region Q(tau): all bins with PDF >= kappa
  central_region_indices <- which(pdf_flat >= kappa)

  # Define the outlier region as the complement of the central region
  outlier_indices <- setdiff(seq_along(pdf_flat), central_region_indices)

  # ---------------------------
  # 4. Prepare data for plotting the outlier bins
  # ---------------------------
  highlight_data <- data.frame(
    x = grid$x[outlier_indices],
    y = grid$y[outlier_indices],
    z = grid$z[outlier_indices],
    value = pdf_flat[outlier_indices]
  )

  # Normalize PDF values for marker sizing
  normalized_pdf <- highlight_data$value / max(highlight_data$value, na.rm = TRUE)

  # ---------------------------
  # 5. Plot the outlier bins using plotly
  # ---------------------------
  fig <- plot_ly(
    data = highlight_data,
    x = ~x,
    y = ~y,
    z = ~z,
    type = "scatter3d",
    mode = "markers",
    marker = list(
      size = ~normalized_pdf * 75,  # Adjust the scaling factor as needed
      color = ~value,
      colorscale = "Viridis",
      showscale = TRUE
    ),
    text = ~paste("PDF Value:", round(value, 4))
  ) %>% layout(
    scene = list(
      xaxis = list(title = "Variable 1 (pr)"),
      yaxis = list(title = "Variable 2 (tas)"),
      zaxis = list(title = "Variable 3 (psl)")
    ),
    title = paste("Outlier Bins (10% of Mass) for Grid Point (Lon:", lon_index, ", Lat:", lat_index, ")")
  )

  # Display the plot
  fig

}

# Assume GC_result11$label_attribution is your label matrix (dimensions: height x width)
label_matrix <- GC_result11$label_attribution

# Number of rows in the label matrix
n_rows <- nrow(label_matrix)

# Initialize a vector to store the mean difference for each adjacent pair of rows
row_diff_means <- numeric(n_rows - 1)

# Loop over each adjacent row pair (comparing row i and row i+1)
for(i in 1:(n_rows - 1)) {
  # diff_vec is TRUE where labels differ between row i and row i+1
  diff_vec <- label_matrix[i, ] != label_matrix[i+1, ]
  # Compute the fraction of pixels that differ
  row_diff_means[i] <- mean(diff_vec)
}

# Now compare the first row with the last row (wrap-around)
wrap_diff <- mean(label_matrix[1, ] != label_matrix[n_rows, ])

# Combine the differences from adjacent rows and the wrap-around
all_diffs <- c(row_diff_means, wrap_diff)

# Plot a histogram of the differences for each row pair (including wrap-around)
hist(all_diffs,
     main = "Histogram of Differences Between Adjacent Rows (with wrap-around)",
     xlab = "Fraction of pixels with differing labels",
     col = "lightblue",
     border = "gray")