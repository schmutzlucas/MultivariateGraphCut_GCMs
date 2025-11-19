# ======================================================================
#  compute_nd_pdf_multi()
#
#  Computes:
#   • 3-D joint PDFs (nbins3d^3)
#   • 2-D pairwise flat PDFs (nbins2d^2)
#   • 1-D marginal PDFs (nbins1d)
#   • per-cell means
#
#  Returns a list with:
#    pdf3$present, pdf3$future
#    pdf2$present, pdf2$future
#    pdf1$present, pdf1$future
#    mean$present, mean$future
# ======================================================================
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

  mean_pres <- array(NA, c(nlon,nlat,nmods))
  mean_fut  <- array(NA, c(nlon,nlat,nmods))

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

    # normalize
    tot_p3 <- rowSums(p3_p, dims=2); p3_p<-sweep(p3_p,1:2,tot_p3,"/")
    tot_f3 <- rowSums(p3_f, dims=2); p3_f<-sweep(p3_f,1:2,tot_f3,"/")

    norm2d <- function(a) sweep(a,1:2,rowSums(a,dims=2),"/")
    for(k in pdf2_names){
      p2_p[[k]]<-norm2d(p2_p[[k]])
      p2_f[[k]]<-norm2d(p2_f[[k]])
    }

    norm1d <- function(a) {
      # a: [lon, lat, nbins1d]
      tot <- rowSums(a, dims = 2)          # sum over bins → [lon, lat]
      a_norm <- sweep(a, 1:2, tot, "/")    # divide each (lon,lat,bin) by tot[lon,lat]
      a_norm[!is.finite(a_norm)] <- 0      # protect against 0/0
      a_norm
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


  # ── return ────────────────────────────────────────────────────────────
  list(
    pdf3 = list(present = pdf3_pres, future = pdf3_fut),
    pdf2 = list(present = pdf2_pres, future = pdf2_fut),
    pdf1 = list(present = pdf1_pres, future = pdf1_fut),
    mean = list(present = mean_pres,  future = mean_fut )
  )

}
