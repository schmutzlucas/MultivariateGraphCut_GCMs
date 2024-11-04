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
#' @param pdf_ref An array representing the kde_ref dataset for the optimization.
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
GraphCutHellinger_xD <- function(
  pdf_models_future,
  h_dist,
  weight_data,
  weight_smooth,
  nBins,
  seed,
  verbose,
  rebuild
){

  n_labs      <- length(model_names)
  width       <- ncol(pdf_models_future[,,,1])
  height      <- nrow(pdf_models_future[,,,1])
  print(width)
  print(height)


  # Permuting longitude and latitude since the indexing isn't the same in R and in C++
  # changed: c(aperm(sum_h_dist, c(2, 1, 3))) call was redundant
  # when go from matrix to vector
  h_dist_cpp <- c(aperm(h_dist, c(2, 1, 3)))
  pdf_models_cpp <- c(aperm(pdf_models_future, c(3, 2, 1, 4)))
  # pdf_models_cpp <- c(pdf_models_future)


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
    rebuild = rebuild, showOutput = FALSE, verbose = FALSE
  )

  # cat("Creating SmoothCost function...  ")
  # ptrSmoothCost <- cppXPtr(
  #   code = 'float smoothFn(int p1, int p2, int l1, int l2, Rcpp::List extraData)
  #   {
  #     int nbVariables = extraData["n_variables"];
  #     int numPix = extraData["numPix"];
  #     float weight = extraData["weight"];
  #     NumericVector data = extraData["data"];
  #     int totPix = numPix;
  #
  #     float cost = 0;
  #
  #     for (int k = 0; k < nbVariables; k++) {
  #       cost += std::abs(data[k + (p1 * nbVariables + totPix * l1)] - data[k + (p2 * nbVariables + totPix * l2)]);
  #     }
  #
  #     return(weight * cost);
  #   }',
  #   includes = c("#include <math.h>", "#include <Rcpp.h>"),
  #   rebuild = TRUE, showOutput = FALSE, verbose = FALSE
  # )

  cat("Creating SmoothCost function...  ")
  ptrSmoothCost <- cppXPtr(
    code = 'float smoothFn(int p1, int p2, int l1, int l2, Rcpp::List extraData)
{
    int numPix = extraData["numPix"];
    float weight = extraData["weight"];
    NumericVector data = extraData["data"];
    int nBins = extraData["nBins"];

    float cost  = 0.0f;
    float tmp1  = 0.0f;
    float tmp2  = 0.0f;
    float diff1 = 0.0f;
    float diff2 = 0.0f;
    int index1, index2;

    // Precompute common terms to reduce redundant calculations
    int offset_p1_l1 = (p1 + numPix * l1) * nBins;
    int offset_p1_l2 = (p1 + numPix * l2) * nBins;
    int offset_p2_l1 = (p2 + numPix * l1) * nBins;
    int offset_p2_l2 = (p2 + numPix * l2) * nBins;

    for (int i = 0; i < nBins; i++) {
        // Compute indices for data access
        index1 = offset_p1_l1 + i;
        index2 = offset_p1_l2 + i;

        // Compute the difference of square roots
        diff1 = sqrt(data[index1]) - sqrt(data[index2]);

        // Accumulate the squared difference
        tmp1 += diff1 * diff1;

        // Repeat for p2
        index1 = offset_p2_l1 + i;
        index2 = offset_p2_l2 + i;

        diff2 = sqrt(data[index1]) - sqrt(data[index2]);
        tmp2 += diff2 * diff2;
    }

    // Compute the Hellinger distances and total cost
    cost = (sqrt(tmp1) + sqrt(tmp2)) / sqrt(2.0f);

    // Debugging: Print each individual cost
    // Rcpp::Rcout << "Cost between (" << p1 << ", " << p2 << ") with labels (" << l1 << ", " << l2 << "): " << weight * cost << std::endl;

    return(weight * cost);
}
',
    includes = c("#include <math.h>", "#include <Rcpp.h>"),
    rebuild = rebuild, showOutput = FALSE, verbose = FALSE
  )


  # Creation of the data and smooth cost
  gco$setDataCost(ptrDataCost, list(numPix  = width * height,
                                    data    = h_dist_cpp,
                                    weight  = weight_data))

  gco$setSmoothCost(ptrSmoothCost, list(numPix  = width * height,
                                        data    = pdf_models_cpp,
                                        weight  = weight_smooth,
                                        nBins   = nBins))


  # # Initialization matrix based on the best model (h_dist)
  # mae_list <- numeric(n_labs)
  # for(i in seq_along(mae_list)){
  #   mae_list[[i]] <- mean(abs(h_dist[,,i]))
  # }
  # best_label <- which.min(mae_list)-1 # in C++ label indices start at 0
  # print(best_label)
  # for(z in 0:((width*height)-1)){
  #   # Label is set as the best average model
  #   gco$setLabel(z, best_label)
  # }

  # Initializing randomly
  set.seed(seed)
  for(z in 0:((width*height)-1)){
    random_label <- sample(0:(n_labs-1), 1) # Sample a random index uniformly
    gco$setLabel(z, random_label)
  }

  # Print the data energy with descriptive text
  cat("Data Energy of Current Labeling:\n")
  cat("  Total data energy: ", gco$giveDataEnergy(), "\n\n")

  # Print the smooth energy with descriptive text
  cat("Smooth Energy of Current Labeling:\n")
  cat("  Total smoothness energy: ", gco$giveSmoothEnergy(), "\n\n")

  # for(z in 0:(length(width*height)-1)){
  #  gco$setLabel(z, 0)
  # }

  # Optimizing the MRF energy with alpha-beta swap
  # -1 refers to the optimization until convergence
  cat("Starting GraphCut optimization...  ")
  print(format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  begin <- Sys.time()
  gco_exec <- gco$swap(-1)
  print(gco_exec)
  time_spent <- Sys.time()-begin
  print(format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  cat("GraphCut optimization done :  ")
  print(time_spent)

  data_cost         <- gco$giveDataEnergy()
  smooth_cost       <- gco$giveSmoothEnergy()
  data_smooth_list  <- list("Data cost" = data_cost, "Smooth cost" = smooth_cost)

  label_attribution <- matrix(0, nrow = height, ncol = width)
  h_dist_GC_present <- matrix(0, nrow = height, ncol = width)
  for(j in 1:height){
    for(i in 1:width){
      label_attribution[j,i] <- gco$whatLabel((i - 1) + width * (j - 1)) ### Permuting from the C++ indexing to the R indexing
    }
  }

  label_attribution <- label_attribution + 1


  tmp <- array(NA, dim = dim(label_attribution))
  for(l in 0:(length(model_names))){
    islabel <- which(label_attribution == l)
    tmp[islabel] <- h_dist[,,(l)][islabel]
  }

  h_dist_GC_present <- tmp


  # gc_result <- vector("list",length=2)
  gc_result <- list("label_attribution" = label_attribution,
                    "Data and smooth cost" = data_smooth_list,
                    'h_dist_GC_present' = h_dist_GC_present)

  return(gc_result)
}