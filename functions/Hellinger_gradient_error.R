gradient_hdist <- function(hdist_map) {
  mask <- !is.na(hdist_map)
  gradient_error <- matrix(0, nrow(hdist_map), ncol(hdist_map))
  denom <- matrix(0, nrow(hdist_map), ncol(hdist_map))

  ileft   <- 1:(nrow(hdist_map) - 1)
  iright  <- 2:nrow(hdist_map)
  itop    <- 1:(ncol(hdist_map) - 1)
  ibottom <- 2:ncol(hdist_map)

  # Only count neighbors where both cells are not NA
  valid_lr <- mask[ileft, ] & mask[iright, ]
  valid_tb <- mask[, itop] & mask[, ibottom]

  # left–right
  diff_lr <- abs(hdist_map[ileft, ] - hdist_map[iright, ])
  diff_lr[!valid_lr] <- 0

  gradient_error[ileft, ][valid_lr]  <- gradient_error[ileft, ][valid_lr]  + diff_lr[valid_lr]
  gradient_error[iright, ][valid_lr] <- gradient_error[iright, ][valid_lr] + diff_lr[valid_lr]

  denom[ileft, ][valid_lr]  <- denom[ileft, ][valid_lr]  + 1
  denom[iright, ][valid_lr] <- denom[iright, ][valid_lr] + 1

  # top–bottom
  diff_tb <- abs(hdist_map[, itop] - hdist_map[, ibottom])
  diff_tb[!valid_tb] <- 0

  gradient_error[, itop][valid_tb]    <- gradient_error[, itop][valid_tb]    + diff_tb[valid_tb]
  gradient_error[, ibottom][valid_tb] <- gradient_error[, ibottom][valid_tb] + diff_tb[valid_tb]

  denom[, itop][valid_tb]    <- denom[, itop][valid_tb]    + 1
  denom[, ibottom][valid_tb] <- denom[, ibottom][valid_tb] + 1

  gradient <- gradient_error / denom
  gradient[denom == 0] <- NA
  gradient
}
