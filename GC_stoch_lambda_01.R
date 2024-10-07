# Stochastic Graph cuts test
# load('202404071152_my_workspace_ERA5_allmodels_stoch_lambda01.RData')

# Install and load necessary libraries
list_of_packages <- read.table("package_list.txt", sep="\n")$V1
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages))
  install.packages(new.packages, repos = "https://cloud.r-project.org")

library(devtools)
lapply(list_of_packages, library, character.only = TRUE)
install_github("thaos/gcoWrapR")


# Loading local functions
source_code_dir <- 'functions/' #The directory where all functions are saved.
file_paths <- list.files(source_code_dir, full.names = T)
for(path in file_paths){source(path)}


N_IT <- 5  # Number of iterations

# Initialize the list to store results outside the loop
GC_result_hellinger_lambda01_5_s <- vector("list", N_IT)
error_list <- vector("list", N_IT)  # To store error messages, if any

lambda <- 0.1

weight_data <- 1 - lambda
weight_smooth <- lambda


for(i in 1:N_IT) {
  # Wrap the function call in tryCatch to handle errors
  GC_result_hellinger_lambda01_5_s[[i]] <- tryCatch({
    # Call the GraphCutHellinger2D_new2 function and store the result
    GraphCutHellinger2D_new2(pdf_ref = kde_ref,
                             kde_models = kde_models,
                             pdf_models_future = kde_models_future,
                             h_dist = h_dist,
                             weight_data = weight_data,
                             weight_smooth = weight_smooth,
                             nBins = nbins1d^2,
                             verbose = TRUE,
                             rebuild = FALSE)
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

# Get the current date and time
current_time <- Sys.time()

# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")

# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "_my_workspace_ERA5_allmodels_stoch_lambda01_5s.RData")

# Save the workspace using the generated filename
save.image(file = filename, compress = FALSE)



avg_h_dist <- list()
data_cost <- c()
smooth_cost <- c()
for (i in seq_along(GC_result_hellinger_lambda01_20)) {
  avg_h_dist[[i]] <- mean(c(GC_result_hellinger_lambda01_20[[i]]$h_dist))
  data_cost[i]    <- GC_result_hellinger_lambda01_20[[i]]$`Data and smooth cost`$`Data cost`
  smooth_cost[i]  <- GC_result_hellinger_lambda01_20[[i]]$`Data and smooth cost`$`Smooth cost`
}


GC_res_lambda01_20_matrix <- array(dim = c(length(lon), length(lat), 20))


for(i in 1:20) {
  GC_res_lambda01_20_matrix[,,i] <- GC_result_hellinger_lambda01_20[[i]]$label_attribution
}

# Initialize the matrix to store entropy values
entropy_matrix <- array(NaN, c(length(lon), length(lat)))

# Loop over each pixel location
for (i in 1:360) {
  for (j in 1:181) {
    # Extract labels for the current pixel across all realizations
    labels <- GC_res_lambda01_20_matrix[i, j, ]

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


# Initialize the matrix to store the count of unique models
unique_model_count_matrix <- array(NA, c(length(lon), length(lat)))

# Loop over each pixel location
for (i in 1:360) {
  for (j in 1:181) {
    # Extract labels for the current pixel across all realizations
    labels <- GC_res_lambda01_20_matrix[i, j, ]

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


h_dist_future_lambda01 <- list()
avg_hdist_future_lambda01 <- list()
for (i in 1:20) {
  tmp <- array(NA, dim = dim(GC_result_hellinger_new[[1]]$label_attribution))
  for(j in seq_along(variables)){
    for(l in 0:(length(model_names))){
      islabel <- which(GC_result_hellinger_lambda01_20[[i]]$label_attribution == l)
      tmp[islabel] <- h_dist_future[,,(l)][islabel]

    }
  }
  h_dist_future_lambda01[[i]] <- tmp
  avg_hdist_future_lambda01[[i]] <- mean(h_dist_future_lambda01[[i]])
}
