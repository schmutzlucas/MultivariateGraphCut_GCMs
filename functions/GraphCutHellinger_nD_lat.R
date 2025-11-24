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
GraphCutHellinger_nD_lat <- function(
  pdf_models_future,
  h_dist,
  weight_data,
  weight_smooth,
  nBins,
  lat,
  seed,
  verbose,
  rebuild
) {
  n_labs <- dim(h_dist)[3]
  height <- ncol(pdf_models_future[,,,1])
  width <- nrow(pdf_models_future[,,,1])

  print(height)
  print(width)

  if (length(lat) != height) {
    stop("Error: Latitude vector length does not match grid height.")
  }

  lat_weights <- matrix(rep(abs(cos(lat * pi / 180)), each = width),
                      nrow = width, ncol = height, byrow = FALSE)

  print(dim(lat_weights))

  # Permuting the arrays for C++ indexing
  h_dist_cpp <- c(aperm(h_dist, c(1, 2, 3)))
  pdf_models_cpp <- c(aperm(pdf_models_future, c(3, 1, 2, 4)))
  lat_weights_cpp <-  c(aperm(lat_weights, c(1, 2)))


  # Instantiate the GraphCut environment
  gco <- new(GCoptimizationGridGraph, width, height, n_labs)

  print(gco)

  # DataCost function using 2D latitude weights
  ptrDataCost <- cppXPtr(
    code = 'float dataFn(int p, int l, Rcpp::List extraData)
{
    int width         = extraData["width"];
    int height        = extraData["height"];
    int numPix        = width * height;
    float weight_global = extraData["weight"];
    Rcpp::NumericVector data        = extraData["data"];
    Rcpp::NumericVector lat_weights = extraData["lat_weights"];

    float lat_weight = lat_weights[p];
    float dval       = data[p + numPix * l];

    // For debugging: only print for p=0, 10000, 20000,... or some special condition
    if (p % 10000 == 0 && l == 3) {
        Rcpp::Rcout << "DEBUG dataFn(): p=" << p
                    << ", l=" << l
                    << ", lat_weight=" << lat_weight
                    << ", dval=" << dval
                    << std::endl;
    }

    // Return the cost
    return(weight_global * dval * lat_weight);
}',
    includes = c("#include <math.h>", "#include <Rcpp.h>"),
    rebuild = rebuild, showOutput = FALSE, verbose = FALSE
  )

  cat("Creating SmoothCost function...  ")
  ptrSmoothCost <- cppXPtr(
    code = 'float smoothFn(int p1, int p2, int l1, int l2, Rcpp::List extraData)
{
    int width = extraData["width"];
    int height = extraData["height"];
    int numPix = width * height;
    float weight_global = extraData["weight"];
    NumericVector data = extraData["data"];
    NumericVector lat_weights = extraData["lat_weights"];
    int nBins = extraData["nBins"];

    // Ensure p1 and p2 are within the valid range
    if (p1 < 0 || p1 >= numPix || p2 < 0 || p2 >= numPix) {
        Rcpp::stop("Pixel index out of bounds in smoothFn.");
    }

    // Retrieve latitude weights safely
    float lat_weight = (lat_weights[p1] + lat_weights[p2]) / 2.0;

    float cost = 0.0f;
    float tmp1 = 0.0f;
    float tmp2 = 0.0f;
    float diff1, diff2;
    int index1, index2;

    int offset_p1_l1 = (p1 + numPix * l1) * nBins;
    int offset_p1_l2 = (p1 + numPix * l2) * nBins;
    int offset_p2_l1 = (p2 + numPix * l1) * nBins;
    int offset_p2_l2 = (p2 + numPix * l2) * nBins;

    for (int i = 0; i < nBins; i++) {
        index1 = offset_p1_l1 + i;
        index2 = offset_p1_l2 + i;
        diff1 = sqrt(data[index1]) - sqrt(data[index2]);
        tmp1 += diff1 * diff1;

        index1 = offset_p2_l1 + i;
        index2 = offset_p2_l2 + i;
        diff2 = sqrt(data[index1]) - sqrt(data[index2]);
        tmp2 += diff2 * diff2;
    }

    cost = (sqrt(tmp1) + sqrt(tmp2)) / sqrt(2.0f);
    return weight_global * lat_weight * cost;
    // return weight_global *  cost;
}
',
    includes = c("#include <math.h>", "#include <Rcpp.h>"),
    rebuild = rebuild, showOutput = TRUE, verbose = FALSE
  )


  # Set DataCost and SmoothCost
  gco$setDataCost(ptrDataCost, list(
    numPix = width * height,
    width = width,
    height = height,
    data = h_dist_cpp,
    weight = weight_data,
    lat_weights = lat_weights_cpp
  ))

  gco$setSmoothCost(ptrSmoothCost, list(
    width = width,
    height = height,
    data = pdf_models_cpp,
    weight = weight_smooth,
    nBins = nBins,
    lat_weights = lat_weights_cpp
  ))


  # Initialize labels randomly
  # set.seed(seed)
  # for (z in 0:((width * height) - 1)) {
  #   random_label <- sample(0:(n_labs - 1), 1)
  #   gco$setLabel(z, random_label)
  # }

  # Initialize labels pseudo-randomly

  RNGkind(kind = "L'Ecuyer-CMRG", normal.kind = "Inversion") # For reproducibility across OS
  set.seed(seed)

  init_labels <- sample.int(n_labs, width * height, replace = TRUE) - 1
  for (z in 0:((width * height) - 1)) {
    gco$setLabel(z, init_labels[z + 1])
  }


  # Perform graph cut optimization
  cat("Starting GraphCut optimization...  ")
  print(format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  begin <- Sys.time()
  gco$swap(-1)  # Run until convergence
  time_spent <- Sys.time() - begin
  print(format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  cat("GraphCut optimization done :  ")
  print(time_spent)

  # Retrieve results
  data_cost <- gco$giveDataEnergy()
  smooth_cost <- gco$giveSmoothEnergy()
  data_smooth_list <- list("Data cost" = data_cost, "Smooth cost" = smooth_cost)

  # Extract label attribution
  label_attribution <- matrix(0, nrow = width, ncol = height)
  for (j in 1:height) {
    for (i in 1:width) {
      label_attribution[i, j] <- gco$whatLabel((i - 1) + width * (j - 1))
    }
  }
  label_attribution <- label_attribution + 1

  # Clean up and force garbage collection
  rm(gco)
  gc()

  # Return results
  list(
    "label_attribution" = label_attribution,
    "Data and smooth cost" = data_smooth_list
  )
}
