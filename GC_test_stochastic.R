# Stochastic Graph cuts tests


N_IT <- 100

# Initialize the list to store results outside the loop
GC_result_hellinger <- vector("list", N_IT)

for(i in 1:N_IT) {
  # Graphcut hellinger labelling and store the result indexed by i
  GC_result_hellinger[[i]] <- GraphCutHellinger2D(kde_ref = kde_ref,
                                                  kde_models = kde_models,
                                                  models_smoothcost = models_matrix_nrm$future,
                                                  weight_data = 1,
                                                  weight_smooth = 1,
                                                  verbose = TRUE)
}