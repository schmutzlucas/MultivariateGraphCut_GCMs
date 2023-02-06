list.of.packages <- c("RcppXPtrUtils","devtools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages,repos = "http://cran.us.r-project.org")
lapply(list.of.packages, library, character.only = TRUE)
# install_github("thaos/gcoWrapR")
library(gcoWrapR)
# TODO documentation of gc
#' @title
#' Graph cut optimization
#'
#' @description
#' This function produces a map of labels where each grid-point is affected to
#' one model. It uses gcoWrapR (https://github.com/thaos/gcoWrapR) and c++ based
#' gco-v3.0 (https://vision.cs.uwaterloo.ca/code/).
#'
#' @param data should be an array
#' @param method allows the user the chose the normalization method.
#' Currently: Standard Score or Min Max
#'
#' @examples
#' Normalize(precipitation, StdSc)
#' Normalize(temperature, MinMax)
#'
#' @return
#' Returns the results of the graphcut as an array of labels
#'
GraphCutOptimization <- function(
  reference,
  models_datacost,
  models_smoothcost,
  weight_data,
  weight_smooth
){

  n_labs      <- length(models_datacost[1, 1, , 1])
  n_variables <- length(reference[1, 1, 1, ])
  width       <- ncol(reference)
  height      <- nrow(reference)

  print(height)
  print(width)
  print(n_labs)
  print(n_variables)


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
    includes = c("#include <math.h>", "#include <Rcpp.h>"),
    rebuild = FALSE, showOutput = FALSE, verbose = FALSE
  )

  cat("Creating SmoothCost function...  ")
  #TODO sortie dans un fichier
  ptrSmoothCost <- cppXPtr(
    code = 'float smoothFn(int p1, int p2, int l1, int l2, Rcpp::List extraData)
{
  int nbVariables        = extraData["n_variables"];
  int numPix             = extraData["numPix"];
  float weight           = extraData["weight"];
  NumericVector data     = extraData["data"];
  int totPix             = numPix * nbVariables;

  float cost = 0;

  for (int k = 0; k < nbVariables; k++) {
    cost += std::abs(data[k + (p1 * nbVariables + totPix * l1)] -
              data[k + (p1 * nbVariables + totPix * l2)]) +
              std::abs(data[k + (p2 * nbVariables + totPix * l1)] -
              data[k + (p2 * nbVariables + totPix * l2)]);
  }

  return(weight * cost);
}',
    includes = c("#include <math.h>", "#include <Rcpp.h>"),
    rebuild = FALSE, showOutput = FALSE, verbose = FALSE
  )

  # Preparing the data to perform GraphCut
  # TODO Test
  bias <- array(0, c(height, width, n_labs))
  for (i in 1:n_labs) {
    for (j in 1:n_variables) {
      bias[,,i] <- bias[,,i] + abs(models_datacost[,, i, j] - reference[[j]])
    }
  }

  # Permuting longitude and latitude since the indexing isn't the same in R and in C++
  # changed: c(aperm(bias, c(2, 1, 3))) call was redundant
  # when go from matrix to vector
  bias_cpp <- c(aperm(bias, c(2, 1, 3)))
  smooth_cpp <- c(aperm(models_smoothcost, c(4, 1,2, 3)))

  # Creation of the data and smooth cost
  gco$setDataCost(ptrDataCost, list(numPix = width * height,
                                    data = bias_cpp,
                                    weight = weight_data))

  gco$setSmoothCost(ptrSmoothCost, list(numPix  = width * height,
                                        data = smooth_cpp,
                                        weight = weight_smooth,
                                        n_variables = n_variables))

  # Creating the initialization matrix based on the best model (bias)
  # TODO Implement random version?
  mae_list <- numeric(n_labs)
  for(i in seq_along(mae_list)){
    mae_list[[i]] <- mean(abs(bias[,,i]))
  }
  best_label <- which.min(mae_list)-1 # in C++ label indices start at 0
  for(z in 0:(length(reference)-1)){
    gco$setLabel(z, best_label)
  }

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

  gc_result <- vector("list",length=2)
  gc_result <- list("label_attribution" = label_attribution, "Data and smooth cost" = data_smooth_list)

  return(gc_result)
}
