gradient_mae <- function(ref, var){
  mae <- matrix(0, nrow(ref), ncol(ref))
  denom <- matrix(0, nrow(ref), ncol(ref))

  ileft <- 1:(nrow(ref)-1)
  iright <- 2:nrow(ref)
  itop <- 1:(ncol(ref)-1)
  ibottom <- 2:ncol(ref)


  denom[ileft, ] <- denom[ileft, ] + 1
  denom[iright, ] <- denom[iright, ] + 1
  denom[, itop] <- denom[, itop] + 1
  denom[, ibottom] <- denom[, ibottom] + 1

  mae[ileft, ] <- mae[ileft,] + abs(var[ileft, ] - var[iright, ] - ref[ileft, ] + ref[iright, ])
  mae[iright, ] <- mae[iright,] + abs(var[iright, ] - var[ileft, ] - ref[iright, ] + ref[ileft, ])
  mae[, itop] <- mae[, itop] + abs(var[, itop] - var[, ibottom] - ref[, itop] + ref[, ibottom ])
  mae[, ibottom] <- mae[, ibottom] + abs(var[, ibottom] - var[, itop] - ref[, ibottom] + ref[, itop])

  return(mae / denom)
}

