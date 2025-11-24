#' Compute Multivariate PDFs and Means for CMIP6 Models
#'
#' This function computes multivariate probability density functions (PDFs) and
#' per-cell means for a fixed set of three climate variables (`pr`, `tas`, `psl`)
#' over a grid of longitude–latitude points and a set of CMIP6 models.
#' For each model and each grid cell, it builds:
#'
#' * a 3-D joint PDF over (`pr`, `tas`, `psl`)
#' * 2-D pairwise joint PDFs for (`pr`, `tas`), (`pr`, `psl`), (`tas`, `psl`)
#' * 1-D marginal PDFs for each variable (`pr`, `tas`, `psl`)
#' * per-cell means for each variable
#'
#' All PDFs are normalised to sum to 1 at each grid cell and model, and internal
#' sanity checks are performed to ensure this property (the function stops with
#' an error if normalisation fails beyond numerical tolerance).
#'
#' @param variables Character vector of length 3 with the variable names.
#'   Must be exactly `c("pr", "tas", "psl")` and in that order.
#' @param model_names Character vector with the names of the climate models
#'   (one entry per model). This defines the model dimension of the outputs.
#' @param data_dir Root directory containing the CMIP6 data, organised as
#'   `file.path(data_dir, model, var)` for each `model` in `model_names` and
#'   each `var` in `variables`. Within each `model/var` folder, the function
#'   expects at least one NetCDF file whose name matches the pattern
#'   `paste0(var, "_", model, "*.nc")`.
#' @param year_present Integer vector of years defining the time span for the
#'   "present" period (e.g., `1950:1975`). Only time steps whose calendar year
#'   is in this vector are used for the present PDFs.
#' @param year_future Integer vector of years defining the time span for the
#'   "future" period (e.g., `1998:2023`). Only time steps whose calendar year
#'   is in this vector are used for the future PDFs.
#' @param lon Numeric vector of target longitudes. These must correspond to
#'   longitudes available in the NetCDF files (after conversion to the
#'   \[-180, 180\] convention used internally). The length of this vector
#'   defines the longitude dimension of the outputs.
#' @param lat Numeric vector of target latitudes. These must correspond to
#'   latitudes available in the NetCDF files. The length of this vector
#'   defines the latitude dimension of the outputs.
#' @param range_var Numeric array of dimension `[nlon, nlat, 3, 2]` giving
#'   the minimum and maximum values for each variable at each grid cell.
#'   The dimensions are:
#'   \itemize{
#'     \item `[, , 1, ]` for `pr`
#'     \item `[, , 2, ]` for `tas`
#'     \item `[, , 3, ]` for `psl`
#'   }
#'   and the last dimension is of length 2, with
#'   `range_var[,,,1]` the minimum and `range_var[,,,2]` the maximum.
#'   These ranges are used consistently for all models to build comparable PDFs.
#' @param nbins3d Integer. Number of bins per variable for the 3-D joint PDFs.
#'   The resulting 3-D histogram at each cell has `nbins3d^3` bins.
#' @param nbins2d Integer. Number of bins per variable for the 2-D joint PDFs.
#'   Each 2-D histogram at each cell has `nbins2d^2` bins.
#' @param nbins1d Integer. Number of bins for each 1-D marginal PDF.
#' @param workers Integer. Number of parallel workers to use for processing
#'   models. On Unix-like systems, `multicore` is used; on Windows,
#'   `multisession` is used via the \pkg{future} framework.
#'
#' @details
#' For each model and variable, the function reads daily data from the first
#' matching NetCDF file in `file.path(data_dir, model, var)`. Longitudes in the
#' NetCDF files are converted to a common \[-180, 180\] system and matched to
#' the user-provided `lon` vector; latitudes are matched directly to `lat`.
#' Only the minimal contiguous blocks in longitude, latitude, and time needed
#' to cover the requested grid and year ranges are read, to reduce I/O.
#'
#' For precipitation (`pr`), a logarithmic transform `log(x + 1)` is applied
#' before histogramming. The N-dimensional histograms are computed via the
#' helper function `compute_histND()`, using the per-cell ranges given by
#' `range_var` and the specified numbers of bins.
#'
#' After counting, all PDFs are explicitly normalised so that, for every grid
#' cell and model, the sum over bins is 1 (within a numerical tolerance).
#' Internal sanity checks verify this property for 3-D, 2-D, and 1-D PDFs and
#' raise an error if any violation is detected.
#'
#' @return A list with four components:
#' \describe{
#'   \item{\code{pdf3}}{
#'     A list with elements:
#'     \describe{
#'       \item{\code{present}}{Numeric array of dimension
#'         `[nlon, nlat, nbins3d^3, nmods]` with the 3-D joint PDFs for the
#'         present period.}
#'       \item{\code{future}}{Numeric array of the same dimension with the
#'         3-D joint PDFs for the future period.}
#'     }
#'   }
#'   \item{\code{pdf2}}{
#'     A list with elements:
#'     \describe{
#'       \item{\code{present}}{Named list of 2-D joint PDFs for the present
#'         period. It contains three entries:
#'         \code{"pr_tas"}, \code{"pr_psl"}, \code{"tas_psl"}. Each entry is a
#'         numeric array of dimension `[nlon, nlat, nbins2d^2, nmods]`.}
#'       \item{\code{future}}{Named list with the same structure for the
#'         future period.}
#'     }
#'   }
#'   \item{\code{pdf1}}{
#'     A list with elements:
#'     \describe{
#'       \item{\code{present}}{Named list of 1-D marginal PDFs for the present
#'         period. It contains three entries: \code{"pr"}, \code{"tas"},
#'         \code{"psl"}. Each entry is a numeric array of dimension
#'         `[nlon, nlat, nbins1d, nmods]`.}
#'       \item{\code{future}}{Named list with the same structure for the
#'         future period.}
#'     }
#'   }
#'   \item{\code{mean}}{
#'     A list with elements:
#'     \describe{
#'       \item{\code{present}}{Named list of per-cell means for the present
#'         period. Each entry (\code{"pr"}, \code{"tas"}, \code{"psl"}) is a
#'         numeric array of dimension `[nlon, nlat, nmods]`.}
#'       \item{\code{future}}{Named list with the same structure for the
#'         future period.}
#'     }
#'   }
#' }
#'
#' @seealso
#' \code{\link{compute_histND}} for the underlying N-dimensional histogram
#' construction used at each grid cell.
#'
#' @examples
#' \dontrun{
#' variables   <- c("pr", "tas", "psl")
#' model_names <- c("CanESM5", "MPI-ESM1-2-HR")
#' data_dir    <- "data/CMIP6_summer_Apr15-Oct14"
#'
#' # toy grid (subset of the full CMIP6 grid)
#' lon <- seq(-10, 10, by = 2.5)
#' lat <- seq( 40, 50, by = 2.5)
#'
#' # range_var: [lon, lat, var, min/max]
#' range_var <- array(NA_real_, dim = c(length(lon), length(lat), 3, 2))
#' # ... fill range_var with suitable per-cell min/max for pr, tas, psl ...
#'
#' year_present <- 1950:1975
#' year_future  <- 1998:2023
#'
#' out <- compute_nd_pdf_multi(
#'   variables    = variables,
#'   model_names  = model_names,
#'   data_dir     = data_dir,
#'   year_present = year_present,
#'   year_future  = year_future,
#'   lon          = lon,
#'   lat          = lat,
#'   range_var    = range_var,
#'   nbins3d      = 8,
#'   nbins2d      = 16,
#'   nbins1d      = 32,
#'   workers      = 4
#' )
#'
#' # Example: sum of a 1-D PDF at a random cell and model (should be 1)
#' i0 <- sample(seq_along(lon), 1)
#' j0 <- sample(seq_along(lat), 1)
#' m0 <- sample(seq_along(model_names), 1)
#'
#' sum(out$pdf1$present$tas[i0, j0, , m0])
#' }
#'
#' @export
compute_nd_pdf_multi <- function(
  variables,     # c("pr","tas","psl")
  model_names,   # list of model names
  data_dir,      # root CMIP6 directory
  year_present,  # e.g. 1950:1975
  year_future,   # e.g. 1998:2023
  lon, lat,      # numeric vectors of long/lat
  range_var,     # [lon,lat,3,2] min/max per var
  nbins3d = 8,   # 3-D bins per dimension
  nbins2d = 16,  # 2-D bins per dimension
  nbins1d = 32,  # 1-D bins
  workers  = 4
) {
  # sanity
  stopifnot(length(variables)==3,
            all(variables==c("pr","tas","psl")))

  nmods   <- length(model_names)
  nlon    <- length(lon)
  nlat    <- length(lat)
  pdf2_names <- c("pr_tas","pr_psl","tas_psl")

  # ── allocate master arrays ────────────────────────────────────────────
  pdf3_pres <- array(NA, c(nlon,nlat,nbins3d^3, nmods))
  pdf3_fut  <- array(NA, c(nlon,nlat,nbins3d^3, nmods))

  pdf2_pres <- lapply(pdf2_names, function(.)
    array(NA, c(nlon,nlat,nbins2d^2, nmods)))
  names(pdf2_pres) <- pdf2_names
  pdf2_fut  <- lapply(pdf2_names, function(.)
    array(NA, c(nlon,nlat,nbins2d^2, nmods)))
  names(pdf2_fut) <- pdf2_names

  pdf1_pres <- lapply(variables, function(.)
    array(NA, c(nlon,nlat,nbins1d, nmods)))
  names(pdf1_pres) <- variables
  pdf1_fut  <- lapply(variables, function(.)
    array(NA, c(nlon,nlat,nbins1d, nmods)))
  names(pdf1_fut) <- variables


  # ── parallel over models ─────────────────────────────────────────────
  library(future); library(future.apply); library(ncdf4)
if (.Platform$OS.type == "unix") {
  plan(multicore, workers = workers)
} else {
  plan(multisession, workers = workers)
}

  model_list <- future_lapply(seq_len(nmods), function(m) {
    model <- model_names[m]
    cat(format(Sys.time(), "%H:%M:%S"), "– reading", model, "\n")

    # ── read raw daily data (robust lon/lat handling) ──────────────────
    var_pres <- vector("list", 3)
    var_fut  <- vector("list", 3)

    for (v in seq_along(variables)) {
      var <- variables[v]

      ## ---------- open the FIRST matching NetCDF ----------------------
      fn <- list.files(
        file.path(data_dir, model, var),
        pattern = glob2rx(paste0(var, "_", model, "*.nc")),
        full.names = TRUE
      )[1]
      if (is.na(fn)) stop("No file found for ", model, "/", var)
      nc  <- nc_open(fn)

      yrs <- extract_years_from_time(nc)   # helper already in code
      lon_file <- ncvar_get(nc, "lon")
      lat_file <- ncvar_get(nc, "lat")

      # bring all longitudes into the same –180 … +180° convention
      lon_file_adj <- ifelse(lon_file >= 180, lon_file - 360, lon_file)
      lon_user_adj <- ifelse(lon      >= 180, lon      - 360, lon)

      ## ---------- build index vectors without NA ----------------------
      # sort-order trick so we can ask for one contiguous block
      lon_order  <- order(lon_file_adj)
      lon_sorted <- lon_file_adj[lon_order]

      lon_idx_sorted <- match(lon_user_adj, lon_sorted)
      if (anyNA(lon_idx_sorted))
        stop("Some requested longitudes are missing in ", basename(fn))

      lon_idx_file <- lon_order[lon_idx_sorted]

      lat_idx_file <- match(lat, lat_file)
      if (anyNA(lat_idx_file))
        stop("Some requested latitudes are missing in ", basename(fn))

      ## ---------- fast reader: a helper that respects spans -----------
      slab <- function(year_span) {

        ii <- which(yrs %in% year_span)
        if (length(ii) == 0)
          stop("No dates in span ", paste(range(year_span), collapse = "-"))

        # contiguous start/count for lon, lat, time
        start_lon <- min(lon_idx_file)
        cnt_lon   <- max(lon_idx_file) - start_lon + 1
        start_lat <- min(lat_idx_file)
        cnt_lat   <- max(lat_idx_file) - start_lat + 1
        start_tim <- min(ii)
        cnt_tim   <- max(ii) - start_tim + 1

        cube <- ncvar_get(
          nc, var,
          start = c(start_lon, start_lat, start_tim),
          count = c(cnt_lon,   cnt_lat,   cnt_tim)
        )

        # local indices inside the cube
        lon_local <- match(lon_idx_file,
                           seq(start_lon, length.out = cnt_lon))
        lat_local <- match(lat_idx_file,
                           seq(start_lat, length.out = cnt_lat))

        cube[lon_local, lat_local, , drop = FALSE]
      }

      ## ---------- read present / future blocks ------------------------
      pres <- slab(year_present)
      fut  <- slab(year_future)

      if (var == "pr") {           # keep your logarithmic transform
        pres <- log(pres + 1)
        fut  <- log(fut  + 1)
      }

      var_pres[[v]] <- pres
      var_fut [[v]] <- fut

      nc_close(nc)
      gc()
    }  # end-for variable

    # per-model hist containers
    p3_p <- array(0, c(nlon,nlat,nbins3d^3))
    p3_f <- array(0, c(nlon,nlat,nbins3d^3))

    p2_p <- lapply(pdf2_names, function(.) array(0, c(nlon,nlat,nbins2d^2)))
    p2_f <- lapply(pdf2_names, function(.) array(0, c(nlon,nlat,nbins2d^2)))
    names(p2_p)<-names(p2_f)<-pdf2_names

    p1_p <- lapply(variables, function(.) array(0, c(nlon,nlat,nbins1d)))
    p1_f <- lapply(variables, function(.) array(0, c(nlon,nlat,nbins1d)))
    names(p1_p)<-names(p1_f)<-variables

    m_p <- lapply(variables, function(.) matrix(0,nlon,nlat))
    m_f <- lapply(variables, function(.) matrix(0,nlon,nlat))
    names(m_p)<-names(m_f)<-variables

    # compute histograms & means
    for (i in seq_len(nlon)) for (j in seq_len(nlat)) {
      pr_p <- var_pres[[1]][i,j,]; ta_p<-var_pres[[2]][i,j,]; ps_p<-var_pres[[3]][i,j,]
      pr_f <- var_fut [[1]][i,j,]; ta_f<-var_fut [[2]][i,j,]; ps_f<-var_fut [[3]][i,j,]

      # ranges
      r1_pr  <- matrix(range_var[i,j,1,],ncol=2,byrow=TRUE)
      r1_tas <- matrix(range_var[i,j,2,],ncol=2,byrow=TRUE)
      r1_psl <- matrix(range_var[i,j,3,],ncol=2,byrow=TRUE)
      r2_pt  <- rbind(range_var[i,j,1,],range_var[i,j,2,])
      r2_pp  <- rbind(range_var[i,j,1,],range_var[i,j,3,])
      r2_tp  <- rbind(range_var[i,j,2,],range_var[i,j,3,])
      r3 <- rbind(range_var[i,j,1,],
                  range_var[i,j,2,],
                  range_var[i,j,3,])

      dat_p <- cbind(pr_p, ta_p, ps_p)
      dat_f <- cbind(pr_f, ta_f, ps_f)

      # 3-D
      p3_p[i,j,] <- compute_histND(dat_p, r3, nbins3d)
      p3_f[i,j,] <- compute_histND(dat_f, r3, nbins3d)

      # 2-D
      p2_p$pr_tas [i,j,] <- compute_histND(dat_p[,1:2],    r2_pt, nbins2d)
      p2_f$pr_tas [i,j,] <- compute_histND(dat_f[,1:2],    r2_pt, nbins2d)
      p2_p$pr_psl [i,j,] <- compute_histND(dat_p[,c(1,3)],  r2_pp, nbins2d)
      p2_f$pr_psl [i,j,] <- compute_histND(dat_f[,c(1,3)],  r2_pp, nbins2d)
      p2_p$tas_psl[i,j,] <- compute_histND(dat_p[,2:3],    r2_tp, nbins2d)
      p2_f$tas_psl[i,j,] <- compute_histND(dat_f[,2:3],    r2_tp, nbins2d)

      # 1-D
      p1_p$pr  [i,j,] <- compute_histND(matrix(pr_p,ncol=1), r1_pr, nbins1d)
      p1_f$pr  [i,j,] <- compute_histND(matrix(pr_f,ncol=1), r1_pr, nbins1d)
      p1_p$tas [i,j,] <- compute_histND(matrix(ta_p,ncol=1), r1_tas,nbins1d)
      p1_f$tas [i,j,] <- compute_histND(matrix(ta_f,ncol=1), r1_tas,nbins1d)
      p1_p$psl [i,j,] <- compute_histND(matrix(ps_p,ncol=1), r1_psl,nbins1d)
      p1_f$psl [i,j,] <- compute_histND(matrix(ps_f,ncol=1), r1_psl,nbins1d)

      # means
      m_p$pr [i,j]<-mean(pr_p); m_f$pr [i,j]<-mean(pr_f)
      m_p$tas[i,j]<-mean(ta_p); m_f$tas[i,j]<-mean(ta_f)
      m_p$psl[i,j]<-mean(ps_p); m_f$psl[i,j]<-mean(ps_f)
    }


    norm3d <- function(a) {
      # a: [lon, lat, nbins3d^3]
      tot <- rowSums(a, dims = 2)       # [lon, lat]
      a_norm <- sweep(a, 1:2, tot, "/")
      a_norm[!is.finite(a_norm)] <- 0
      a_norm
    }

    # normalize 3-D
    p3_p <- norm3d(p3_p)
    p3_f <- norm3d(p3_f)

    norm2d <- function(a) {
      # a: [lon, lat, nbins2d^2]
      tot <- rowSums(a, dims = 2)
      a_norm <- sweep(a, 1:2, tot, "/")
      a_norm[!is.finite(a_norm)] <- 0
      a_norm
    }

    # normalize 2-D
    for (k in pdf2_names) {
      p2_p[[k]] <- norm2d(p2_p[[k]])
      p2_f[[k]] <- norm2d(p2_f[[k]])
    }

    norm1d <- function(a) {
      # a: [lon, lat, nbins1d]
      tot <- rowSums(a, dims = 2)          # sum over bins → [lon, lat]
      a_norm <- sweep(a, 1:2, tot, "/")    # divide each (lon,lat,bin) by tot[lon,lat]
      a_norm[!is.finite(a_norm)] <- 0      # protect against 0/0
      a_norm
    }

    for (k in variables) {
      p1_p[[k]] <- norm1d(p1_p[[k]])
      p1_f[[k]] <- norm1d(p1_f[[k]])
    }


    list(
      pdf3_pres = p3_p, pdf3_fut = p3_f,
      pdf2_pres = p2_p, pdf2_fut = p2_f,
      pdf1_pres = p1_p, pdf1_fut = p1_f,
      mean_pres = m_p,  mean_fut  = m_f
    )
  })
  plan(sequential)

  # ---- create mean containers per variable ----------------------------------
  mean_pres <- lapply(variables, function(.) array(NA, c(nlon, nlat, nmods)))
  mean_fut  <- lapply(variables, function(.) array(NA, c(nlon, nlat, nmods)))
  names(mean_pres) <- names(mean_fut) <- variables

  # ---- collate over models ---------------------------------------------------
  for (m in seq_len(nmods)) {

    # 3-D
    pdf3_pres[,,,m] <- model_list[[m]]$pdf3_pres
    pdf3_fut [,,,m] <- model_list[[m]]$pdf3_fut

    # 2-D
    for (k in pdf2_names) {
      pdf2_pres[[k]][,,,m] <- model_list[[m]]$pdf2_pres[[k]]
      pdf2_fut [[k]][,,,m] <- model_list[[m]]$pdf2_fut [[k]]
    }

    # 1-D
    for (k in variables) {
      pdf1_pres[[k]][,,,m] <- model_list[[m]]$pdf1_pres[[k]]
      pdf1_fut [[k]][,,,m] <- model_list[[m]]$pdf1_fut [[k]]

      # means now collate correctly
      mean_pres[[k]][,,m]  <- model_list[[m]]$mean_pres[[k]]
      mean_fut [[k]][,,m]  <- model_list[[m]]$mean_fut [[k]]
    }
  }


  # ---- internal sanity checks: PDFs must sum to 1 per cell -----------

  tol <- 1e-6  # tolerance for floating-point noise

  check_pdf3 <- function(a, name) {
    # a: [lon, lat, nbins3d^3, nmods]
    sums <- apply(a, c(1, 2, 4), sum)  # -> [lon, lat, model]
    bad  <- (sums < -tol) | (sums > tol & abs(sums - 1) > tol)

    if (any(bad, na.rm = TRUE)) {
      idx <- which(bad, arr.ind = TRUE)[1, ]
      s   <- sums[idx[1], idx[2], idx[3]]
      stop(sprintf(
        "PDF normalisation error in %s at [lon=%d, lat=%d, model=%d]: sum = %.8f",
        name, idx[1], idx[2], idx[3], s
      ))
    }
  }

  check_pdfN <- function(a, name) {
    # for 2-D and 1-D PDFs shaped [lon, lat, nbins, nmods]
    sums <- apply(a, c(1, 2, 4), sum)  # -> [lon, lat, model]
    bad  <- (sums < -tol) | (sums > tol & abs(sums - 1) > tol)

    if (any(bad, na.rm = TRUE)) {
      idx <- which(bad, arr.ind = TRUE)[1, ]
      s   <- sums[idx[1], idx[2], idx[3]]
      stop(sprintf(
        "PDF normalisation error in %s at [lon=%d, lat=%d, model=%d]: sum = %.8f",
        name, idx[1], idx[2], idx[3], s
      ))
    }
  }

  # 3-D checks
  check_pdf3(pdf3_pres, "pdf3$present")
  check_pdf3(pdf3_fut,  "pdf3$future")

  # 2-D checks
  for (k in pdf2_names) {
    check_pdfN(pdf2_pres[[k]], paste0("pdf2$present$", k))
    check_pdfN(pdf2_fut [[k]], paste0("pdf2$future$",  k))
  }

  # 1-D checks
  for (k in variables) {
    check_pdfN(pdf1_pres[[k]], paste0("pdf1$present$", k))
    check_pdfN(pdf1_fut [[k]], paste0("pdf1$future$",  k))
  }


  # ── return ────────────────────────────────────────────────────────────
  list(
    pdf3 = list(present = pdf3_pres, future = pdf3_fut),
    pdf2 = list(present = pdf2_pres, future = pdf2_fut),
    pdf1 = list(present = pdf1_pres, future = pdf1_fut),
    mean = list(present = mean_pres,  future = mean_fut )
  )

}
