list.of.packages <- c("RcppXPtrUtils","devtools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages,repos = "http://cran.us.r-project.org")
lapply(list.of.packages, library, character.only = TRUE)
# install_github("thaos/gcoWrapR")
library(gcoWrapR)
#' Graph cut optimization
#'
#' This function performs graph cut optimization using the gco-v3.0 C++ library and gcoWrapR package.
#' It produces a map of labels where each grid-point is assigned to one model.
#'
#' @param kde_ref An array representing the kde_ref dataset for the optimization.
#' @param models_datacost An array representing the models' data cost for the optimization.
#' @param models_smoothcost An array representing the models' smooth cost for the optimization.
#' @param weight_data A numeric value representing the weight for the data cost.
#' @param weight_smooth A numeric value representing the weight for the smooth cost.
#' @param verbose A logical value indicating whether or not to print information during the optimization process.
#'
#' @return A list containing the label attribution matrix, data and smooth cost, and execution time.
#'
#' @examples
#' # Load example data
#' data("example_data")
#'
#' # Perform graph cut optimization
#' GC_result <- GraphCutOptimization(
#' kde_ref = example_data$kde_ref,
#' models_datacost = example_data$models_datacost,
#' models_smoothcost = example_data$models_smoothcost,
#' weight_data = 1,
#' weight_smooth = 1,
#' verbose = TRUE
#' )
#'
#' # Print results
#' print(GC_result)
#'
#' @references
#' https://github.com/thaos/gcoWrapR
#' https://vision.cs.uwaterloo.ca/code/
#'
#' @import gcoWrapR
#' @export
GraphCutHellinger2D_new <- function(
  kde_ref,
  kde_models,
  models_smoothcost,
  weight_data,
  weight_smooth,
  verbose
){

  n_labs      <- length(model_names)
  width       <- ncol(kde_ref)
  height      <- nrow(kde_ref)
  n_variables <- 2



  # Computing the sum of hellinger distances between models and reference --> used as datacost
  h_dist <- array(data = 0, dim = c(length(lon), length(lat),
                                    length(model_names)))

  # Loop through variables and models
  m <- 1
  for (model_name in model_names) {
    for (i in seq_along(lon)) {
      for (j in seq_along(lat)) {
        # Compute Hellinger distance
        h_dist_unchecked <- sqrt(sum((sqrt(kde_models[i, j, , m]) - sqrt(kde_ref[i, j, ]))^2)) / sqrt(2)

        # Replace NaN with 0
        h_dist[i, j, m] <- replace(h_dist_unchecked, is.nan(h_dist_unchecked), 0)
      }
    }
    m <- m + 1
  }


  # Permuting longitude and latitude since the indexing isn't the same in R and in C++
  # changed: c(aperm(sum_h_dist, c(2, 1, 3))) call was redundant
  # when go from matrix to vector
  h_dist_cpp <- c(aperm(h_dist, c(2, 1, 3)))
  smooth_cpp <- c(aperm(models_smoothcost, c(4, 2, 1, 3)))


  # Instanciation of the GraphCut environment

  gco <- new(GCoptimizationGridGraph, width, height, n_labs)

  # Preparing the DataCost and SmoothCost functions of the GraphCut in C++

  cat("Creating DataCost function...  ")

  ptrDataCost <- cppXPtr(
    code = 'float dataFn(int p, int l, Rcpp::List extraData)
    {

      int numPix          = extraData["numPix"];
      float weight        = extraData["weight"];
      NumericVector data  = extraData["data"];

      return(weight * data[p + numPix * l]);
    }',
    includes = c("#include <math.h>", "#include <Rcpp.h>", "#include <iostream>"),
    rebuild = TRUE, showOutput = FALSE, verbose = FALSE
  )

  cat("Creating SmoothCost function...  ")
  ptrSmoothCost <- cppXPtr(
    code = 'float smoothFn(int p1, int p2, int l1, int l2, Rcpp::List extraData) {
  int nbVariables = extraData["n_variables"];
  int numPix = extraData["numPix"];
  float weight = extraData["weight"];
  NumericVector data = extraData["data"];
  int totPix = numPix * nbVariables;

  float cost = 0;

  for (int k = 0; k < nbVariables; k++) {
    cost += std::abs(data[k + (p1 * nbVariables + totPix * l1)] - data[k + (p2 * nbVariables + totPix * l2)]);
  }

  return(weight * cost);
}
',
    includes = c("#include <math.h>", "#include <Rcpp.h>"),
    rebuild = TRUE, showOutput = FALSE, verbose = FALSE
  )




  # Creation of the data and smooth cost
  gco$setDataCost(ptrDataCost, list(numPix = width * height,
                                    data = h_dist_cpp,
                                    weight = weight_data))

  gco$setSmoothCost(ptrSmoothCost, list(numPix  = width * height,
                                        data = smooth_cpp,
                                        weight = weight_smooth,
                                        n_variables = n_variables))

  # Creating the initialization matrix based on the best model (sum_h_dist)
  # # TODO Implement random version?
  mae_list <- numeric(n_labs)
  for(i in seq_along(mae_list)){
    mae_list[[i]] <- mean(abs(h_dist[,,i]))
  }
  best_label <- which.min(mae_list)-1 # in C++ label indices start at 0
  print(best_label)
  for(z in 0:((width*height)-1)){
    # Label is set as the best average model
    # gco$setLabel(z, best_label)

    random_label <- sample(0:(n_labs-1), 1) # Sample a random index uniformly
    gco$setLabel(z, random_label)
    # #   # gco$setLabel(z, -1)
    # #   gco$setLabel(z, 1)
  }

  # for(z in 0:(length(width*height)-1)){
  #  gco$setLabel(z, 7)
  # }

  # Optimizing the MRF energy with alpha-beta swap
  # -1 refers to the optimization until convergence
  cat("Starting GraphCut optimization...  ")
  begin <- Sys.time()
  gco$swap(-1)
  time_spent <- Sys.time()-begin
  cat("GraphCut optimization done :  ")
  print(time_spent)

  data_cost         <- gco$giveDataEnergy()
  smooth_cost       <- gco$giveSmoothEnergy()
  data_smooth_list  <- list("Data cost" = data_cost, "Smooth cost" = smooth_cost)

  label_attribution <- matrix(0,nrow = height,ncol = width)
  for(j in 1:height){
    for(i in 1:width){
      label_attribution[j,i] <- gco$whatLabel((i - 1) + width * (j - 1)) ### Permuting from the C++ indexing to the R indexing
    }
  }

  label_attribution <- label_attribution + 1

  # gc_result <- vector("list",length=2)
  gc_result <- list("label_attribution" = label_attribution, "Data and smooth cost" = data_smooth_list, 'h_dist' = h_dist)

  return(gc_result)
}


