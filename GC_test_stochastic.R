# Stochastic Graph cuts tests

N_IT <- 100  # Number of iterations

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
