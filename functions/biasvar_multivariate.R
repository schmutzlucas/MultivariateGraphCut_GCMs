bias_var <- function(var_future, ref_future_nrm, var_future_nrm, labeling){
  width <- ncol(ref_future_nrm[[1]])
  height <- nrow(ref_future_nrm[[1]])
  nlabs <- length(var_future[[1]][1,1,])
  bias_gc <- list()
  var_gc_nrm <- var_gc <- list(matrix(0,nrow = height,ncol = width),matrix(0,nrow = height,ncol = width))
  nvar <- length(var_future)
  nref <- length(ref_future_nrm)

  for(j in 1:nvar){
    for(l in 0:(nlabs-1)){
      islabel <- which(labeling == l)
      var_gc[[j]][islabel] <- var_future[[j]][,,(l+1)][islabel]
    }
  }
  for(j in 1:nref){
    for(l in 0:(nlabs-1)){
      islabel <- which(labeling == l)
      var_gc_nrm[[j]][islabel] <- var_future_nrm[[j]][,,(l+1)][islabel]
    }
  }
  if(is.null(ref_future_nrm))
  {
    bias_gc <-  NULL
  } else {
      for(j in 1:nref){
        var_gc_nrm[[j]][islabel] <- var_future_nrm[[j]][,,(l+1)][islabel]
        bias_gc[[j]] <- (var_gc_nrm[[j]] -  ref_future_nrm[[j]])
      }
  }
  bias_var <- list("var" = var_gc, "bias" = bias_gc)
  return(bias_var)
}
