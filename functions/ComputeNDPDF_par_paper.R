# ======================================================================
#  compute_nd_pdf_multi()
#
#  - 3-D joint PDF (nbins3d³)
#  - 2-D pairwise PDFs (flat nbins2d²)
#  - 1-D marginal PDFs (nbins1d)
#  - cell‐wise means
# ======================================================================
compute_nd_pdf_multi <- function(variables, model_names, data_dir,
                                 year_present, year_future,
                                 lon, lat, range_var,
                                 nbins3d = 8, nbins2d = 16, nbins1d = 32,
                                 workers  = 4)
{
  stopifnot(length(variables) == 3,
            all(variables == c("pr","tas","psl")))

  nmodels <- length(model_names)
  nlon    <- length(lon)
  nlat    <- length(lat)

  # ── allocate master arrays ────────────────────────────────────────────
  pdf3_pres <- array(NA, c(nlon, nlat, nbins3d^3, nmodels))
  pdf3_fut  <- array(NA, c(nlon, nlat, nbins3d^3, nmodels))

  pdf2_names <- c("pr_tas","pr_psl","tas_psl")
  pdf2_pres  <- lapply(pdf2_names, function(x)
    array(NA, c(nlon, nlat, nbins2d^2, nmodels)))
  pdf2_fut   <- lapply(pdf2_names, function(x)
    array(NA, c(nlon, nlat, nbins2d^2, nmodels)))

  pdf1_pres  <- lapply(variables, function(x)
    array(NA, c(nlon, nlat, nbins1d, nmodels)))
  pdf1_fut   <- lapply(variables, function(x)
    array(NA, c(nlon, nlat, nbins1d, nmodels)))
  names(pdf1_pres) <- names(pdf1_fut) <- variables

  mean_pres  <- lapply(variables, function(x)
    array(NA, c(nlon, nlat, nmodels)))
  mean_fut   <- lapply(variables, function(x)
    array(NA, c(nlon, nlat, nmodels)))
  names(mean_pres) <- names(mean_fut) <- variables

  # ── parallel over models ─────────────────────────────────────────────
  plan(multisession, workers = workers)
  model_list <- future_lapply(seq_along(model_names), function(m) {
    model <- model_names[m]
    cat(format(Sys.time(), "%H:%M:%S"), "– reading", model, "\n")

    # read variables
    var_dat_pres <- vector("list", 3)
    var_dat_fut  <- vector("list", 3)
    for (v in seq_along(variables)) {
      var  <- variables[v]
      fdir <- file.path(data_dir, model, var)
      fn   <- list.files(fdir, pattern = glob2rx(paste0(var, "_", model, "*.nc")),
                         full.names = TRUE)[1]
      nc   <- nc_open(fn)

      yyyy    <- extract_years_from_time(nc)
      lon_idx <- match(lon, ncvar_get(nc, "lon"))
      lat_idx <- match(lat, ncvar_get(nc, "lat"))

      slice <- function(year_span) {
        idx <- which(yyyy %in% year_span)
        ncvar_get(nc, var,
                  start = c(min(lon_idx), min(lat_idx), min(idx)),
                  count = c(length(lon_idx), length(lat_idx), length(idx)))
      }

      pres <- slice(year_present)
      fut  <- slice(year_future)
      if (var == "pr") {
        pres <- log(pres + 1)
        fut  <- log(fut  + 1)
      }

      var_dat_pres[[v]] <- pres
      var_dat_fut [[v]] <- fut
      nc_close(nc)
    }

    # per-model containers
    p3_pres <- array(0, c(nlon, nlat, nbins3d^3))
    p3_fut  <- array(0, c(nlon, nlat, nbins3d^3))

    p2_pres <- lapply(pdf2_names, function(x)
      array(0, c(nlon, nlat, nbins2d^2)))
    p2_fut  <- lapply(pdf2_names, function(x)
      array(0, c(nlon, nlat, nbins2d^2)))
    names(p2_pres) <- names(p2_fut) <- pdf2_names

    p1_pres <- lapply(variables, function(x)
      array(0, c(nlon, nlat, nbins1d)))
    p1_fut  <- lapply(variables, function(x)
      array(0, c(nlon, nlat, nbins1d)))
    names(p1_pres) <- names(p1_fut) <- variables

    m_pres  <- lapply(variables, function(x) matrix(0, nlon, nlat))
    m_fut   <- lapply(variables, function(x) matrix(0, nlon, nlat))
    names(m_pres)  <- names(m_fut) <- variables

    # loop over grid
    for (i in seq_len(nlon)) for (j in seq_len(nlat)) {
      pr_p  <- var_dat_pres[[1]][i,j,];   pr_f <- var_dat_fut[[1]][i,j,]
      tas_p <- var_dat_pres[[2]][i,j,]; tas_f <- var_dat_fut[[2]][i,j,]
      psl_p <- var_dat_pres[[3]][i,j,]; psl_f <- var_dat_fut[[3]][i,j,]

      rng1_pr  <- matrix(range_var[i,j,1,], ncol=2, byrow=TRUE)
      rng1_tas <- matrix(range_var[i,j,2,], ncol=2, byrow=TRUE)
      rng1_psl <- matrix(range_var[i,j,3,], ncol=2, byrow=TRUE)

      rng2_pr_tas  <- rbind(range_var[i,j,1,], range_var[i,j,2,])
      rng2_pr_psl  <- rbind(range_var[i,j,1,], range_var[i,j,3,])
      rng2_tas_psl <- rbind(range_var[i,j,2,], range_var[i,j,3,])

      rng3 <- rbind(range_var[i,j,1,],
                    range_var[i,j,2,],
                    range_var[i,j,3,])

      dat_p <- cbind(pr_p, tas_p, psl_p)
      dat_f <- cbind(pr_f, tas_f, psl_f)

      # 3-D
      p3_pres[i,j,] <- compute_histND(dat_p, rng3, nbins3d)
      p3_fut [i,j,] <- compute_histND(dat_f, rng3, nbins3d)

      # 2-D
      p2_pres$pr_tas [i,j,] <- compute_histND(dat_p[,1:2],    rng2_pr_tas,  nbins2d)
      p2_fut$pr_tas  [i,j,] <- compute_histND(dat_f[,1:2],    rng2_pr_tas,  nbins2d)

      p2_pres$pr_psl [i,j,] <- compute_histND(dat_p[,c(1,3)],  rng2_pr_psl,  nbins2d)
      p2_fut$pr_psl  [i,j,] <- compute_histND(dat_f[,c(1,3)],  rng2_pr_psl,  nbins2d)

      p2_pres$tas_psl[i,j,] <- compute_histND(dat_p[,2:3],    rng2_tas_psl, nbins2d)
      p2_fut$tas_psl [i,j,] <- compute_histND(dat_f[,2:3],    rng2_tas_psl, nbins2d)

      # 1-D
      p1_pres$pr  [i,j,] <- compute_histND(matrix(pr_p, ncol=1),  rng1_pr,  nbins1d)
      p1_fut$pr   [i,j,] <- compute_histND(matrix(pr_f, ncol=1),  rng1_pr,  nbins1d)

      p1_pres$tas [i,j,] <- compute_histND(matrix(tas_p,ncol=1),   rng1_tas, nbins1d)
      p1_fut$tas  [i,j,] <- compute_histND(matrix(tas_f,ncol=1),   rng1_tas, nbins1d)

      p1_pres$psl [i,j,] <- compute_histND(matrix(psl_p,ncol=1),   rng1_psl, nbins1d)
      p1_fut$psl  [i,j,] <- compute_histND(matrix(psl_f,ncol=1),   rng1_psl, nbins1d)

      # means
      m_pres$pr [i,j] <- mean(pr_p);  m_fut$pr [i,j] <- mean(pr_f)
      m_pres$tas[i,j] <- mean(tas_p); m_fut$tas[i,j] <- mean(tas_f)
      m_pres$psl[i,j] <- mean(psl_p); m_fut$psl[i,j] <- mean(psl_f)
    }

    # normalise
    tot_p <- rowSums(p3_pres, dims = 2)
    tot_f <- rowSums(p3_fut,  dims = 2)
    p3_pres <- sweep(p3_pres, 1:2, tot_p, "/")
    p3_fut  <- sweep(p3_fut,  1:2, tot_f, "/")

    norm2d <- function(a) sweep(a,1:2,rowSums(a,dims=2),"/")
    for (k in pdf2_names) {
      p2_pres[[k]] <- norm2d(p2_pres[[k]])
      p2_fut [[k]] <- norm2d(p2_fut [[k]])
    }

    norm1d <- function(a) sweep(a,1:2,rowSums(a),"/")
    for (k in variables) {
      p1_pres[[k]] <- norm1d(p1_pres[[k]])
      p1_fut [[k]] <- norm1d(p1_fut [[k]])
    }

    list(pdf3_pres = p3_pres, pdf3_fut = p3_fut,
         pdf2_pres = p2_pres, pdf2_fut = p2_fut,
         pdf1_pres = p1_pres, pdf1_fut = p1_fut,
         mean_pres = m_pres,  mean_fut = m_fut)
  })  # future_lapply
  plan(sequential)

  # ---------- collate over models ---------------------------------------
  for (m in seq_along(model_names)) {
    # 3-D joint
    pdf3_pres[ , , , m] <- model_list[[m]]$pdf3_pres
    pdf3_fut [ , , , m] <- model_list[[m]]$pdf3_fut

    # 2-D pairs
    for (k in pdf2_names) {
      pdf2_pres[[k]][ , , , m] <- model_list[[m]]$pdf2_pres[[k]]
      pdf2_fut [[k]][ , , , m] <- model_list[[m]]$pdf2_fut [[k]]
    }

    # 1-D marginals
    for (k in variables) {
      pdf1_pres[[k]][ , , , m] <- model_list[[m]]$pdf1_pres[[k]]
      pdf1_fut [[k]][ , , , m] <- model_list[[m]]$pdf1_fut [[k]]

      # and the means are 3-D ([lon, lat, model])
      mean_pres[[k]][ , , m] <- model_list[[m]]$mean_pres[[k]]
      mean_fut [[k]][ , , m] <- model_list[[m]]$mean_fut [[k]]
    }
  }



  list(
    pdf3 = list(present = pdf3_pres, future = pdf3_fut),
    pdf2 = list(present = pdf2_pres, future = pdf2_fut),
    pdf1 = list(present = pdf1_pres, future = pdf1_fut),
    mean = list(present = mean_pres, future = mean_fut)
  )
}
# ======================================================================
