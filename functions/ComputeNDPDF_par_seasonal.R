compute_nd_pdf_multi_seasonal <- function(
  variables,
  model_names,
  data_dir,
  year_present,
  year_future,
  lon, lat,
  range_var,
  season_start_md = c(4, 15),
  season_end_md   = c(10, 14),
  nbins3d = 8,
  nbins2d = 16,
  nbins1d = 32,
  workers = 4
) {
  stopifnot(length(variables) == 3, all(variables == c("pr", "tas", "psl")))
  stopifnot(length(season_start_md) == 2, length(season_end_md) == 2)

  nmods <- length(model_names)
  nlon <- length(lon)
  nlat <- length(lat)
  pdf2_names <- c("pr_tas", "pr_psl", "tas_psl")

  in_annual_window <- function(month, day, start_md, end_md) {
    md <- month * 100 + day
    start_val <- start_md[1] * 100 + start_md[2]
    end_val <- end_md[1] * 100 + end_md[2]

    if (start_val <= end_val) {
      md >= start_val & md <= end_val
    } else {
      md >= start_val | md <= end_val
    }
  }

  extract_time_ymd <- function(nc_obj) {
    time_dim <- if ("time" %in% names(nc_obj$dim)) {
      "time"
    } else if ("valid_time" %in% names(nc_obj$dim)) {
      "valid_time"
    } else {
      stop("No time dimension ('time' or 'valid_time') found.")
    }

    time_raw <- ncdf4::ncvar_get(nc_obj, time_dim)
    time_units <- ncdf4::ncatt_get(nc_obj, time_dim, "units")$value

    if (!is.character(time_units) || !grepl("since", time_units, ignore.case = TRUE)) {
      stop("Unsupported or missing time units in NetCDF: ", time_dim)
    }

    origin_str <- sub(".*since\\s+", "", time_units, ignore.case = TRUE)
    origin_str <- sub(" UTC$", "", origin_str, ignore.case = TRUE)

    origin_time <- suppressWarnings(as.POSIXct(origin_str, tz = "UTC"))
    if (is.na(origin_time)) {
      origin_time <- as.POSIXct(as.Date(origin_str), tz = "UTC")
    }
    if (is.na(origin_time)) {
      stop("Could not parse NetCDF time origin: ", origin_str)
    }

    units_l <- tolower(time_units)
    time_sec <- if (grepl("days", units_l)) {
      time_raw * 86400
    } else if (grepl("hours", units_l)) {
      time_raw * 3600
    } else if (grepl("seconds", units_l)) {
      time_raw
    } else {
      stop("Unsupported time unit (expected days/hours/seconds): ", time_units)
    }

    dt <- origin_time + time_sec
    list(
      year = as.integer(format(dt, "%Y")),
      month = as.integer(format(dt, "%m")),
      day = as.integer(format(dt, "%d"))
    )
  }

  pdf3_pres <- array(NA_real_, c(nlon, nlat, nbins3d^3, nmods))
  pdf3_fut <- array(NA_real_, c(nlon, nlat, nbins3d^3, nmods))

  pdf2_pres <- lapply(pdf2_names, function(.) array(NA_real_, c(nlon, nlat, nbins2d^2, nmods)))
  names(pdf2_pres) <- pdf2_names
  pdf2_fut <- lapply(pdf2_names, function(.) array(NA_real_, c(nlon, nlat, nbins2d^2, nmods)))
  names(pdf2_fut) <- pdf2_names

  pdf1_pres <- lapply(variables, function(.) array(NA_real_, c(nlon, nlat, nbins1d, nmods)))
  names(pdf1_pres) <- variables
  pdf1_fut <- lapply(variables, function(.) array(NA_real_, c(nlon, nlat, nbins1d, nmods)))
  names(pdf1_fut) <- variables

  library(future)
  library(future.apply)
  library(ncdf4)

  if (.Platform$OS.type == "unix") {
    plan(multicore, workers = workers)
  } else {
    plan(multisession, workers = workers)
  }

  model_list <- future_lapply(seq_len(nmods), function(m) {
    model <- model_names[m]
    cat(format(Sys.time(), "%H:%M:%S"), "- reading", model, "\n")

    var_pres <- vector("list", 3)
    var_fut <- vector("list", 3)

    for (v in seq_along(variables)) {
      var <- variables[v]

      fn <- list.files(
        file.path(data_dir, model, var),
        pattern = glob2rx(paste0(var, "_", model, "*.nc")),
        full.names = TRUE
      )[1]
      if (is.na(fn)) stop("No file found for ", model, "/", var)
      nc <- nc_open(fn)

      ymd <- extract_time_ymd(nc)
      season_mask <- in_annual_window(ymd$month, ymd$day, season_start_md, season_end_md)

      lon_file <- ncvar_get(nc, "lon")
      lat_file <- ncvar_get(nc, "lat")

      lon_file_adj <- ifelse(lon_file >= 180, lon_file - 360, lon_file)
      lon_user_adj <- ifelse(lon >= 180, lon - 360, lon)

      lon_order <- order(lon_file_adj)
      lon_sorted <- lon_file_adj[lon_order]

      lon_idx_sorted <- match(lon_user_adj, lon_sorted)
      if (anyNA(lon_idx_sorted)) stop("Some requested longitudes are missing in ", basename(fn))
      lon_idx_file <- lon_order[lon_idx_sorted]

      lat_idx_file <- match(lat, lat_file)
      if (anyNA(lat_idx_file)) stop("Some requested latitudes are missing in ", basename(fn))

      slab <- function(year_span) {
        ii <- which(ymd$year %in% year_span & season_mask)
        if (length(ii) == 0) {
          stop(
            "No dates found for year span ",
            paste(range(year_span), collapse = "-"),
            " and seasonal window ",
            sprintf("%02d-%02d", season_start_md[1], season_start_md[2]),
            " to ",
            sprintf("%02d-%02d", season_end_md[1], season_end_md[2])
          )
        }

        start_lon <- min(lon_idx_file)
        cnt_lon <- max(lon_idx_file) - start_lon + 1
        start_lat <- min(lat_idx_file)
        cnt_lat <- max(lat_idx_file) - start_lat + 1
        start_tim <- min(ii)
        cnt_tim <- max(ii) - start_tim + 1

        cube <- ncvar_get(
          nc, var,
          start = c(start_lon, start_lat, start_tim),
          count = c(cnt_lon, cnt_lat, cnt_tim)
        )

        lon_local <- match(lon_idx_file, seq(start_lon, length.out = cnt_lon))
        lat_local <- match(lat_idx_file, seq(start_lat, length.out = cnt_lat))
        tim_local <- ii - start_tim + 1

        cube[lon_local, lat_local, tim_local, drop = FALSE]
      }

      pres <- slab(year_present)
      fut <- slab(year_future)

      if (var == "pr") {
        pres <- log(pres + 1)
        fut <- log(fut + 1)
      }

      var_pres[[v]] <- pres
      var_fut[[v]] <- fut

      nc_close(nc)
      gc()
    }

    p3_p <- array(0, c(nlon, nlat, nbins3d^3))
    p3_f <- array(0, c(nlon, nlat, nbins3d^3))

    p2_p <- lapply(pdf2_names, function(.) array(0, c(nlon, nlat, nbins2d^2)))
    p2_f <- lapply(pdf2_names, function(.) array(0, c(nlon, nlat, nbins2d^2)))
    names(p2_p) <- names(p2_f) <- pdf2_names

    p1_p <- lapply(variables, function(.) array(0, c(nlon, nlat, nbins1d)))
    p1_f <- lapply(variables, function(.) array(0, c(nlon, nlat, nbins1d)))
    names(p1_p) <- names(p1_f) <- variables

    m_p <- lapply(variables, function(.) matrix(0, nlon, nlat))
    m_f <- lapply(variables, function(.) matrix(0, nlon, nlat))
    names(m_p) <- names(m_f) <- variables

    for (i in seq_len(nlon)) for (j in seq_len(nlat)) {
      pr_p <- var_pres[[1]][i, j, ]
      ta_p <- var_pres[[2]][i, j, ]
      ps_p <- var_pres[[3]][i, j, ]
      pr_f <- var_fut[[1]][i, j, ]
      ta_f <- var_fut[[2]][i, j, ]
      ps_f <- var_fut[[3]][i, j, ]

      r1_pr <- matrix(range_var[i, j, 1, ], ncol = 2, byrow = TRUE)
      r1_tas <- matrix(range_var[i, j, 2, ], ncol = 2, byrow = TRUE)
      r1_psl <- matrix(range_var[i, j, 3, ], ncol = 2, byrow = TRUE)
      r2_pt <- rbind(range_var[i, j, 1, ], range_var[i, j, 2, ])
      r2_pp <- rbind(range_var[i, j, 1, ], range_var[i, j, 3, ])
      r2_tp <- rbind(range_var[i, j, 2, ], range_var[i, j, 3, ])
      r3 <- rbind(range_var[i, j, 1, ], range_var[i, j, 2, ], range_var[i, j, 3, ])

      dat_p <- cbind(pr_p, ta_p, ps_p)
      dat_f <- cbind(pr_f, ta_f, ps_f)

      p3_p[i, j, ] <- compute_histND(dat_p, r3, nbins3d)
      p3_f[i, j, ] <- compute_histND(dat_f, r3, nbins3d)

      p2_p$pr_tas[i, j, ] <- compute_histND(dat_p[, 1:2], r2_pt, nbins2d)
      p2_f$pr_tas[i, j, ] <- compute_histND(dat_f[, 1:2], r2_pt, nbins2d)
      p2_p$pr_psl[i, j, ] <- compute_histND(dat_p[, c(1, 3)], r2_pp, nbins2d)
      p2_f$pr_psl[i, j, ] <- compute_histND(dat_f[, c(1, 3)], r2_pp, nbins2d)
      p2_p$tas_psl[i, j, ] <- compute_histND(dat_p[, 2:3], r2_tp, nbins2d)
      p2_f$tas_psl[i, j, ] <- compute_histND(dat_f[, 2:3], r2_tp, nbins2d)

      p1_p$pr[i, j, ] <- compute_histND(matrix(pr_p, ncol = 1), r1_pr, nbins1d)
      p1_f$pr[i, j, ] <- compute_histND(matrix(pr_f, ncol = 1), r1_pr, nbins1d)
      p1_p$tas[i, j, ] <- compute_histND(matrix(ta_p, ncol = 1), r1_tas, nbins1d)
      p1_f$tas[i, j, ] <- compute_histND(matrix(ta_f, ncol = 1), r1_tas, nbins1d)
      p1_p$psl[i, j, ] <- compute_histND(matrix(ps_p, ncol = 1), r1_psl, nbins1d)
      p1_f$psl[i, j, ] <- compute_histND(matrix(ps_f, ncol = 1), r1_psl, nbins1d)

      m_p$pr[i, j] <- mean(pr_p)
      m_f$pr[i, j] <- mean(pr_f)
      m_p$tas[i, j] <- mean(ta_p)
      m_f$tas[i, j] <- mean(ta_f)
      m_p$psl[i, j] <- mean(ps_p)
      m_f$psl[i, j] <- mean(ps_f)
    }

    norm3d <- function(a) {
      tot <- rowSums(a, dims = 2)
      a_norm <- sweep(a, 1:2, tot, "/")
      a_norm[!is.finite(a_norm)] <- 0
      a_norm
    }
    p3_p <- norm3d(p3_p)
    p3_f <- norm3d(p3_f)

    norm2d <- function(a) {
      tot <- rowSums(a, dims = 2)
      a_norm <- sweep(a, 1:2, tot, "/")
      a_norm[!is.finite(a_norm)] <- 0
      a_norm
    }
    for (k in pdf2_names) {
      p2_p[[k]] <- norm2d(p2_p[[k]])
      p2_f[[k]] <- norm2d(p2_f[[k]])
    }

    norm1d <- function(a) {
      tot <- rowSums(a, dims = 2)
      a_norm <- sweep(a, 1:2, tot, "/")
      a_norm[!is.finite(a_norm)] <- 0
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
      mean_pres = m_p, mean_fut = m_f
    )
  })

  plan(sequential)

  mean_pres <- lapply(variables, function(.) array(NA_real_, c(nlon, nlat, nmods)))
  mean_fut <- lapply(variables, function(.) array(NA_real_, c(nlon, nlat, nmods)))
  names(mean_pres) <- names(mean_fut) <- variables

  for (m in seq_len(nmods)) {
    pdf3_pres[, , , m] <- model_list[[m]]$pdf3_pres
    pdf3_fut[, , , m] <- model_list[[m]]$pdf3_fut

    for (k in pdf2_names) {
      pdf2_pres[[k]][, , , m] <- model_list[[m]]$pdf2_pres[[k]]
      pdf2_fut[[k]][, , , m] <- model_list[[m]]$pdf2_fut[[k]]
    }

    for (k in variables) {
      pdf1_pres[[k]][, , , m] <- model_list[[m]]$pdf1_pres[[k]]
      pdf1_fut[[k]][, , , m] <- model_list[[m]]$pdf1_fut[[k]]
      mean_pres[[k]][, , m] <- model_list[[m]]$mean_pres[[k]]
      mean_fut[[k]][, , m] <- model_list[[m]]$mean_fut[[k]]
    }
  }

  list(
    pdf3 = list(present = pdf3_pres, future = pdf3_fut),
    pdf2 = list(present = pdf2_pres, future = pdf2_fut),
    pdf1 = list(present = pdf1_pres, future = pdf1_fut),
    mean = list(present = mean_pres, future = mean_fut)
  )
}
