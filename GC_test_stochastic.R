# Stochastic Graph cuts tests


# Loading local functions
source_code_dir <- 'functions/' #The directory where all functions are saved.
file_paths <- list.files(source_code_dir, full.names = T)
for(path in file_paths){source(path)}

N_IT <- 5  # Number of iterations

# Initialize the list to store results outside the loop
GC_result_hellinger <- vector("list", N_IT)
error_list <- vector("list", N_IT)  # To store error messages, if any

for(i in 1:N_IT) {
  # Wrap the function call in tryCatch to handle errors
  GC_result_hellinger[[i]] <- tryCatch({
    # Attempt to call the GraphCutHellinger2D function
    GraphCutHellinger2D(kde_ref = kde_ref,
                        kde_models = kde_models,
                        models_smoothcost = models_matrix_nrm$future,
                        weight_data = 1,
                        weight_smooth = 1,
                        verbose = TRUE)
  }, error = function(e) {
    # If an error occurs, save the error message and return NULL for this iteration
    error_list[[i]] <- paste("Error in iteration", i, ":", e$message)
    NULL  # Returning NULL to indicate failure for this iteration
  })

  # Optionally, log the error message
  if (!is.null(error_list[[i]])) {
    cat(error_list[[i]], "\n")  # Print the error message to the console
  }
}

# After the loop, you can check which iterations failed
failed_iterations <- which(sapply(error_list, function(x) !is.null(x)))
if(length(failed_iterations) > 0) {
  cat("Iterations that failed:", paste(failed_iterations, collapse=", "), "\n")
} else {
  cat("All iterations completed successfully.\n")
}


for(i in 1:25) {
  hist(GC_result_hellinger[[i]]$label_attribution, nclass = length(model_names))
}




p <- list()
for(i in 1:N_IT) {
  GC_labels <- GC_result_hellinger[[i]]$label_attribution

  label_df <-melt(GC_labels, c("lon", "lat"),
                  value.name = "label_attribution")
  label_df$lat <- label_df$lat -91

  p[[i]] <- ggplot() +
    geom_tile(data=label_df, aes(x=lon, y=lat, fill= factor(label_attribution)),)+
    scale_fill_hue(na.value=NA,
                   labels = model_names,
                   name = 'Model')+
    ggtitle('Label attribution for GC hybrid')+
    # guide = guide_colourbar(
    #   barwidth = 0.4,
    #   ticks.colour = 'black',
    #   ticks.linewidth = 1,
    #   frame.colour = 'black',
    #   draw.ulim = TRUE,
    #   draw.llim = TRUE),)

    borders("world2", colour = 'black', lwd = 0.12) +
    scale_x_continuous(, expand = c(0, 0)) +
    scale_y_continuous(, expand = c(0,0))+
    theme(legend.position = 'bottom')+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(panel.background = element_blank())+
    xlab('Longitude')+
    ylab('Latitude') +
    labs(fill='[mm/day]')+
    theme_bw()+
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.7, 'cm'), #change legend key height
          legend.key.width = unit(0.2, 'cm'), #change legend key width
          legend.title = element_text(size=9), #change legend title font size
          legend.text = element_text(size=7))+ #change legend text font size
    theme(plot.title = element_text(size=16),
          axis.text=element_text(size=7),
          axis.title=element_text(size=9),)+
    easy_center_title()
  p[[i]]
}

p[[1]]
p[[2]]
p[[3]]
p[[4]]
p[[5]]


########################################################################################################################


# Can it be parrallelized?
library(foreach)
library(doParallel)

N_IT <- 10  # Number of iterations

# Register parallel backend to use multiple cores
no_cores <- detectCores() - 1  # Reserve one core for system processes
registerDoParallel(no_cores)

# Prepare to store results and errors
results <- vector("list", N_IT)
errors <- vector("list", N_IT)

# Parallel loop using foreach
results <- foreach(i = 1:N_IT, .errorhandling = "pass") %dopar% {
  tryCatch({
    # Your function call
    GraphCutHellinger2D(kde_ref = kde_ref,
                        kde_models = kde_models,
                        models_smoothcost = models_matrix_nrm$future,
                        weight_data = 1,
                        weight_smooth = 1,
                        verbose = TRUE)
  }, error = function(e) {
    # Error handling
    list(error = paste("Error in iteration", i, ":", e$message))
  })
}

# Optionally process results and errors after the loop
for (i in seq_along(results)) {
  if (is.list(results[[i]]) && !is.null(results[[i]]$error)) {
    errors[[i]] <- results[[i]]$error
    results[[i]] <- NULL
  }
}

# Check errors
if (length(errors[!sapply(errors, is.null)]) > 0) {
  cat("Errors occurred in the following iterations:\n")
  print(errors[!sapply(errors, is.null)])
} else {
  cat("All iterations completed without errors.\n")
}

# Stop the parallel backend
stopImplicitCluster()





########################################################################################################################
# Loading local functions
source_code_dir <- 'functions/' #The directory where all functions are saved.
file_paths <- list.files(source_code_dir, full.names = T)
for(path in file_paths){source(path)}


N_IT <- 1

GC_result_hellinger_new <- GraphCutHellinger2D_stoch_new(kde_ref = kde_ref,
                                                 kde_models = kde_models,
                                                 models_smoothcost = models_matrix_nrm$future,
                                                 weight_data = 1,
                                                 weight_smooth = 1,
                                                 N_IT = N_IT,
                                                 verbose = TRUE)

GC_result_matrix_new <- array(dim = c(length(lon), length(lat), N_IT))

for(i in 1:N_IT) {
  GC_result_matrix_new[,,i] <- GC_result_hellinger_new[[i]]$label_attribution
}


# Get the current date and time
current_time <- Sys.time()

# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")

# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "_my_workspace_ERA5_stochastic_500.RData")

# Save the workspace using the generated filename
save.image(file = filename, compress = FALSE)



# Calculate the maximum value for the x-axis limit
max_value <- max(GC_result_matrix_1000)

# Create a sequence of breaks at every half unit starting from -0.5
breaks <- seq(-0.5, max_value + 0.5, by = 1)

# Generate the histogram with specified breaks
hist(GC_result_matrix_1000, breaks = breaks, xlim = c(0, max_value), xaxt = 'n', freq = FALSE)

# Add custom x-axis ticks centered on the bars
axis(1, at = seq(0, max_value, by = 1), labels = seq(0, max_value, by = 1))


name <- paste0('figure/H_dist_future_GCHellinger_26models_new')
ggsave(paste0(name, '.pdf'), plot = p1, width = 35, height = 25, units = "cm", dpi = 300)
ggsave(paste0(name, '.png'), plot = p1, width = 35, height = 25, units = "cm", dpi = 300)



# Initialize the matrix to store entropy values
entropy_matrix <- array(NaN, c(length(lon), length(lat)))

# Loop over each pixel location
for (i in 1:360) {
 for (j in 1:181) {
   # Extract labels for the current pixel across all realizations
   labels <- GC_result_matrix_1000[i, j, ]

   # Calculate the frequency distribution of the labels
   label_freq <- table(labels) / length(labels)

   # Compute entropy
   entropy <- -sum(label_freq * log2(label_freq))

   # Assign the computed entropy to the entropy matrix
   entropy_matrix[i, j] <- entropy
 }
}


test_df <- melt(entropy_matrix, c("lon", "lat"), value.name = "Entropy")

p1 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Entropy))+
  ggtitle(paste0("Entropy of 1000 GC runs | ", 'Avg', " = ", round(mean(entropy_matrix), 2)))+
  scale_fill_viridis_c(option = "plasma", limits = c(0, max(entropy_matrix, na.rm = TRUE)))+
  borders("world2", colour = 'black', lwd = 0.12) +
  scale_x_continuous(, expand = c(0, 0)) +
  scale_y_continuous(, expand = c(0,0))+
  theme(legend.position = 'bottom')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())+
  xlab('Longitude')+
  ylab('Latitude') +
  labs(fill = expression(H[2]))+
  theme_bw()+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=12))+ #change legend text font size
  theme(plot.title = element_text(size=24),
        plot.subtitle = element_text(size = 20,hjust=0.5),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),)+
  easy_center_title()
p1

name <- paste0('figure/entropy_map_1000')
ggsave(paste0(name, '.pdf'), plot = p1, width = 35, height = 25, units = "cm", dpi = 300)
ggsave(paste0(name, '.png'), plot = p1, width = 35, height = 25, units = "cm", dpi = 300)





# Initialize the matrix to store the count of unique models
unique_model_count_matrix <- array(NA, c(length(lon), length(lat)))

# Loop over each pixel location
for (i in 1:360) {
  for (j in 1:181) {
    # Extract labels for the current pixel across all realizations
    labels <- GC_result_matrix_1000[i, j, ]

    # Count the number of unique labels (models)
    unique_models <- length(unique(labels))

    # Assign the count of unique models to the matrix
    unique_model_count_matrix[i, j] <- unique_models
  }
}

# Calculate the maximum value for the x-axis limit
max_value <- max(unique_model_count_matrix)

# Create a sequence of breaks at every half unit starting from -0.5
breaks <- seq(-0.5, max_value + 0.5, by = 1)

# Generate the histogram with specified breaks
hist(unique_model_count_matrix, breaks = breaks, xlim = c(0, max_value), xaxt = 'n', freq = FALSE, main = "Number of models per grid point", xlab = "", ylab = "")

# Add custom x-axis ticks centered on the bars
axis(1, at = seq(0, max_value, by = 1), labels = seq(0, max_value, by = 1))


min(unique_model_count_matrix)
max(unique_model_count_matrix)



unique_model_count_df <- melt(unique_model_count_matrix, c("lon", "lat"), value.name = "Count")

# Get a color palette with colorRampPalette function
color_palette <- colorRampPalette(c("white", "blue", "#0072B2", "darkblue", "black"))

# Determine the number of unique counts to decide on the number of colors to use
num_colors <- length(unique(unique_model_count_df$Count))

p1 <- ggplot() +
  geom_tile(data=unique_model_count_df, aes(x=lon, y=lat-90, fill=Count))+
  ggtitle(paste0('Per grid-point model count'))+
  scale_fill_viridis_c(option = "plasma", limits = c(0, max(unique_model_count_matrix, na.rm = TRUE)))+
  # guide = guide_colourbar(
  #   barwidth = 0.4,
  #   ticks.colour = 'black',
  #   ticks.linewidth = 1,
  #   frame.colour = 'black',
  #   draw.ulim = TRUE,
  #   draw.llim = TRUE),)
  borders("world2", colour = 'black', lwd = 0.12) +
  scale_x_continuous(, expand = c(0, 0)) +
  scale_y_continuous(, expand = c(0,0))+
  theme(legend.position = 'bottom')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())+
  xlab('Longitude')+
  ylab('Latitude') +
  labs(fill='#M')+
  theme_bw()+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=12))+ #change legend text font size
  theme(plot.title = element_text(size=24),
        plot.subtitle = element_text(size = 20,hjust=0.5),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),)+
  easy_center_title()
p1

name <- paste0('figure/unique_model_count_map_1000')
ggsave(paste0(name, '.pdf'), plot = p1, width = 35, height = 25, units = "cm", dpi = 300)
ggsave(paste0(name, '.png'), plot = p1, width = 35, height = 25, units = "cm", dpi = 300)


# Choose one of the following lines and replace it with the # Insert one of the scale_fill lines here in the above code
# For Viridis palette:
# scale_fill_viridis_c(option = "C", limits = c(0, max(unique_model_count_matrix, na.rm = TRUE)))
# For Spectral palette:
# scale_fill_gradientn(colors = RColorBrewer::brewer.pal(11, "Spectral"), limits = c(0, max(unique_model_count_matrix, na.rm = TRUE)))
# For Plasma palette:
# scale_fill_viridis_c(option = "plasma", limits = c(0, max(unique_model_count_matrix, na.rm = TRUE)))








