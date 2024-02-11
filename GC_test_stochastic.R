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


for(i in 1:N_IT) {
  hist(GC_result_hellinger[[i]]$label_attribution)
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

