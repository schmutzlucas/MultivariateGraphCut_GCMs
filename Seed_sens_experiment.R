# Install and load necessary libraries
list_of_packages <- read.table("package_list.txt", sep="\n")$V1
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages))
  install.packages(new.packages, repos = "https://cloud.r-project.org")

library(devtools)
lapply(list_of_packages, library, character.only = TRUE)
install_github("schmutzlucas/gcoWrapR")


library(future)
plan(sequential)  # no R process parallelism, just native threads



# Loading local functions
source_code_dir <- 'functions/' #The directory where all functions are saved.
file_paths <- list.files(source_code_dir, full.names = T)
for(path in file_paths){source(path)}

load("workspaces/202511221241_workspace_multiRes_3v_gc_results.RData")

# Setting global variables
lon <- -180:179
lat <- -90:90
lon_size <- length(lon)
lat_size <- length(lat)
# Temporal ranges
year_present <<- 1950:1975
year_future <<- 1999:2024
# data directory
data_dir <<- 'data/CMIP6_summer_Apr15-Oct14'


# List of the variable used
variables <- c('pr', 'tas', 'psl')

## 1.  bin-resolution choices
nbins3d <- 8     # 3-D joint
nbins2d <- 16    # 2-D pairs
nbins1d <- 32    # 1-D marginals


# 3) Parameters -----------------------------------------------
smooth_cost <- 0.1
seeds       <- 1:50           # vector of seeds to explore
n_workers   <- 6              # adjust to your machine

# 4) Parallel plan (Linux: multicore) -------------------------
plan(multicore, workers = n_workers)

# 5) One-seed worker function --------------------------------
run_one_seed <- function(this_seed) {

  cat("Running GraphCut for seed =", this_seed, "\n")

  GC_result <- tryCatch({
    GraphCutHellinger_nD_lat(
      pdf_models_future = pdf3_models_fut,
      h_dist            = h_dist_pres,
      weight_data       = 1,
      weight_smooth     = smooth_cost,
      nBins             = nbins_total3d,
      lat               = lat,
      seed              = this_seed,
      verbose           = FALSE,
      rebuild           = TRUE
    )
  }, error = function(e) {
    cat("  ⚠ GraphCut failed for seed", this_seed, ":", e$message, "\n")
    return(NULL)
  })

  if (is.null(GC_result)) {
    return(list(
      seed           = this_seed,
      mean_hdist_fut = NA_real_
    ))
  }

  # -----------------------------------------------------------
  # Compute GC_hdist_pres / GC_hdist_fut for THIS seed
  # (your exact code)
  # -----------------------------------------------------------
  GC_hdist_pres <- matrix(NA_real_, nrow = length(lon), ncol = length(lat))
  GC_hdist_fut  <- matrix(NA_real_, nrow = length(lon), ncol = length(lat))

  for (l in seq_along(model_names)) {  # ensure alignment with model_names
    islabel <- which(GC_result$label_attribution == l)
    if (length(islabel) == 0) next

    GC_hdist_pres[islabel] <- h_dist_pres[ , , l][islabel]
    GC_hdist_fut [islabel] <- h_dist_fut [ , , l][islabel]
  }

  mean_hdist_fut <- mean(GC_hdist_fut, na.rm = TRUE)

  # Optional: manual GC in worker
  gc()

  # Only return compact stats (to keep memory low)
  list(
    seed           = this_seed,
    mean_hdist_fut = mean_hdist_fut
  )
}

# 6) Parallel execution over seeds ----------------------------
# future.seed = TRUE to get reproducible RNG per element
res_list <- future_lapply(seeds, run_one_seed, future.seed = TRUE)

# 7) Assemble results table -----------------------------------
GC_seed_stats <- do.call(
  rbind,
  lapply(res_list, function(x) as.data.frame(x, stringsAsFactors = FALSE))
)

# 8) Save results ---------------------------------------------
out_file <- sprintf("GC_seed_stats_smooth%.3f.rds", smooth_cost)
saveRDS(GC_seed_stats, file = out_file)

cat("\nSaved seed statistics to:", out_file, "\n")

# Optional: quick summary
print(GC_seed_stats)

best_idx  <- which.min(GC_seed_stats$mean_hdist_fut)
best_seed <- GC_seed_stats$seed[best_idx]
cat("\nBest seed:", best_seed,
    "with mean H-dist future =", GC_seed_stats$mean_hdist_fut[best_idx], "\n")


