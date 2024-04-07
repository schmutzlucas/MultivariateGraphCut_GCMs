# load('202404031541_my_workspace_ERA5_allmodels_new.RData')

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


# Hyperparameters
lambdas <- c(0, 0.01, 0.1, 0.3, 0.5, 0.7)

# Precompute h_dist for all models
cat("Starting Hellinger distance computation...  ")
begin <- Sys.time()

# Check if h_dist already exists in the global environment
if (!exists("h_dist")) {
  cat("h_dist does not exist. Starting Hellinger distance computation...  ")

  # Pre-compute the square roots if applicable
  sqrt_kde_models <- sqrt(kde_models)
  sqrt_kde_ref <- sqrt(kde_ref)

  # Initialize h_dist array
  h_dist <- array(data = 0, dim = c(length(lon), length(lat), length(model_names)))

  # Loop to compute Hellinger distance for each cell across all models
  for (m in seq_along(model_names)) {
    for (i in seq_along(lon)) {
      for (j in seq_along(lat)) {
        # Define a function to compute Hellinger distance for a given cell across all models
        compute_hellinger <- function(m, i, j) {
          sqrt(sum((sqrt_kde_models[i, j, , m] - sqrt_kde_ref[i, j, ])^2)) / sqrt(2)
        }

        # Compute Hellinger distance and store in h_dist
        h_dist[i, j, m] <- compute_hellinger(m, i, j)
      }
    }
  }

  # Replace NaN values in h_dist
  h_dist[is.nan(h_dist)] <- 0

  cat("Hellinger distance computation completed in", Sys.time() - begin, "\n")
} else {
  cat("h_dist already exists. Skipping computation.\n")
}
time_spent <- Sys.time()-begin
cat("Hellinger distance computation done in :  ")
print(time_spent)

# Initialize an empty list to store the results
GC_result_hellinger_new <- list()

for (i in seq_along(lambdas)) {
  # Compute the weights based on the current lambdas
  weight_data <- 1 - lambdas[i]
  weight_smooth <- lambdas[i]

  # Use tryCatch to handle errors
  GC_result_hellinger_new[[i]] <- tryCatch({
    # Call the GraphCutHellinger2D_new2 function and store the result
    GraphCutHellinger2D_new2(kde_ref = kde_ref,
                             kde_models = kde_models,
                             kde_models_future = kde_models_future,
                             h_dist = h_dist,
                             weight_data = weight_data,
                             weight_smooth = weight_smooth,
                             nBins = nbins1d^2,
                             verbose = TRUE,
                             rebuild = FALSE)
  }, error = function(e) {
    # Handle the error, e.g., by printing a message and returning NULL or a specific error value
    cat("An error occurred in iteration", i, "with lambda =", lambdas[i], "\nError message:", e$message, "\n")
    # Optionally, return a specific value indicating the error
    NULL  # Or any other value that makes sense in your context
  })

  # Optionally, name each element of the list by its corresponding lambda value for easier reference
  names(GC_result_hellinger_new)[i] <- as.character(lambdas[i])
}

# Get the current date and time
current_time <- Sys.time()

# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")

# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "_my_workspace_ERA5_allmodels_lambda_test.RData")

# Save the workspace using the generated filename
save.image(file = filename, compress = FALSE)


avg_h_dist <- list()
data_cost <- c()
smooth_cost <- c()
for (i in seq_along(lambdas)) {
  avg_h_dist[[i]] <- mean(c(GC_result_hellinger_new[[i]]$h_dist))
  data_cost[i]    <- GC_result_hellinger_new[[i]]$`Data and smooth cost`$`Data cost`
  smooth_cost[i]  <- GC_result_hellinger_new[[i]]$`Data and smooth cost`$`Smooth cost`
}

# Generate the polychrome color palette with 26 colors
color_palette <- pals::glasbey(26)
h <- list()
for (i in seq_along(lambdas)) {
  GC_labels <- GC_result_hellinger_new[[i]]$label_attribution

  label_df <- melt(GC_labels, c("lon", "lat"), value.name = "label_attribution")
  label_df$lat <- label_df$lat - 91

  h[[i]] <- ggplot() +
    geom_tile(data = label_df, aes(x = lon, y = lat, fill = factor(label_attribution))) +
    scale_fill_manual(values = as.vector(color_palette), na.value = NA, guide = FALSE) +  # guide = FALSE to remove legend
    ggtitle('Label attribution for GC hybrid') +
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
          legend.title = element_text(size=16), #change legend title font size
          legend.text = element_text(size=12))+ #change legend text font size
    theme(plot.title = element_text(size=24),
          plot.subtitle = element_text(size = 20,hjust=0.5),
          axis.text=element_text(size=14),
          axis.title=element_text(size=16),)+
    easy_center_title()

}

h_dist_future_lambdas <- list()
avg_hdist_future <- list()
for (i in seq_along(lambdas)) {
  tmp <- array(NA, dim = dim(GC_result_hellinger_new[[1]]$label_attribution))
  for(j in seq_along(variables)){
    for(l in 0:(length(model_names))){
      islabel <- which(GC_result_hellinger_new[[i]]$label_attribution == l)
      tmp[islabel] <- h_dist_future[,,(l)][islabel]

    }
  }
  h_dist_future_lambdas[[i]] <- tmp
  avg_hdist_future[[i]] <- mean(h_dist_future_lambdas[[i]])
}

GC_hellinger_projections_new <- list()

for (i in seq_along(lambdas)) {
  tmp4 <- list()
  j <- 1
  for(var in variables){
    for(l in 0:(length(model_names))){
      islabel <- which(GC_result_hellinger_new[[i]]$label_attribution == l)
      tmp4[[var]][islabel] <- models_matrix$future[,,l,j][islabel]
    }
    j <- j + 1
  }
  GC_hellinger_projections_new[[i]] <- list("tas" = matrix(tmp4$tas, nrow = 360), "pr" = matrix(tmp4$pr * 86400, nrow = 360))
}




tmp_array <- array(data = NA, dim = c(length(lon), length(lat), nbins1d^2))

for(j in 1:length(lon)){
  for(i in 1:length(lat)){
    label <- GC_result_hellinger_new$`0.1`$label_attribution[j,i]
    tmp_array[j, i, ] <- kde_models_future[j, i, , label]
  }
}

GC_hellinger_projections_pdf_01 <- tmp_array




h_dist_map <- array(NA, dim = dim(GC_result_hellinger_new[[3]]$label_attribution))

for(j in seq_along(variables)){
  for(l in 0:(length(model_names))){
    islabel <- which(GC_result_hellinger_new[[3]]$label_attribution == l)
    h_dist_map[islabel] <- h_dist_future[,,(l)][islabel]
  }
}


test_df <- melt(h_dist_map, c("lon", "lat"), value.name = "Bias")

p5 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias))+
labs(subtitle = 'Projection period : 2000 - 2022')+
  ggtitle(paste0('GC Hellinger', ': Mean Hellinger distance = ', round(mean(h_dist_map), 2)))+
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

name <- paste0('figure/H_dist_future_GC_hellinger_lambda01_26models')
ggsave(paste0(name, '.pdf'), plot = p5, width = 35, height = 25, units = "cm", dpi = 300)
ggsave(paste0(name, '.png'), plot = p5, width = 35, height = 25, units = "cm", dpi = 300)
