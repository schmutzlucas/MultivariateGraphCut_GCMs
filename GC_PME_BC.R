# Install and load necessary libraries
list_of_packages <- read.table("package_list.txt", sep="\n")$V1
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages))
  install.packages(new.packages, repos = "https://cloud.r-project.org")

library(devtools)
lapply(list_of_packages, library, character.only = TRUE)
# install_github("schmutzlucas/gcoWrapR")

# Loading local functions
source_code_dir <- 'functions/'  # The directory where all functions are saved.
file_paths <- list.files(source_code_dir, full.names = TRUE)
for(path in file_paths){
  source(path)
}


# # Setting global variables
# lon <- -10:120
# lat <- 0:75
# lon_size <- length(lon)
# lat_size <- length(lat)
#
# # Temporal ranges
# year_present <<- 1950:1975
# year_future <<- 2075:2100



# Setting global variables
lon <- -180:179
lat <- -90:90
lon_size <- length(lon)
lat_size <- length(lat)

# Temporal ranges
year_present <<- 1950:1975
year_future <<- 2075:2100

workers <- 3

# Data directory
data_dir <<- 'data/CMIP6_merged_all/'

# List of variables used
variables <- c('pr', 'tas', 'psl')

# Bins for the PDFs (nbins1d is the number of bins per variable)
nbins1d <<- 8
nbins <<- nbins1d^(length(variables))
# (The joint PDF will have nbins1d^n_vars bins)

# Obtain the list of models from a file
model_names <- read.table('model_names_pr_tas_psl_perfect_model.txt')
model_names <- as.list(model_names[['V1']])


# Initialization for multiple GraphCut smooth costs
GC_result_list <- list("0.05" = list(), "0.1" = list(), "0.6" = list(), "1" = list(), "2" = list())
GC_present_list <- GC_result_list
GC_future_list <- GC_result_list
GC_hdist_present_list <- GC_result_list
GC_hdist_future_list <- GC_result_list
GC_partial_hdist_future_list <- GC_result_list

MMM_hdist_present_list <- list()
MMM_hdist_future_list <- list()
MMM_partial_hdist_future_list <- list()
MMM_present_list <- list()
MMM_future_list <- list()

smooth_costs <- c(0.05, 0.1, 0.6, 1, 2)


# Charger les checkpoints si disponibles
checkpoint_files <- list.files("checkpoints", pattern = "\\.rds$", full.names = TRUE)

for (f in checkpoint_files) {
  checkpoint <- readRDS(f)
  list2env(checkpoint, envir = .GlobalEnv)
}



for (m in seq_along(model_names)) {
  model_names <- read.table('model_names_pr_tas_psl_perfect_model.txt')
  model_names <- as.list(model_names[['V1']])

  reference_name <- model_names[[m]]
  model_names <- model_names[-m]

  # Sauter si déjà traité
  if (reference_name %in% names(MMM_hdist_future_list)) {
    cat("Reference", reference_name, "already processed. Skipping.\n")
    next
  }

  cat("Processing model", reference_name, "as reference\n")

  time_optimized <- system.time({
    results <- compute_nd_pdf_bias_corrected_2(
      variables, reference_name, model_names, data_dir,
      year_present, year_future, lon, lat, nbins1d,
      workers = workers, buffer = 0.10, verbose = TRUE
    )
  })
  cat("Time taken for compute_nd_pdf_bias_corrected:",
      format_time(time_optimized["elapsed"]), "\n")

  pdf_ref_present <- results$pdf_ref$present
  pdf_models_present <- results$pdf_models$present
  pdf_ref_future <- results$pdf_ref$future
  pdf_models_future <- results$pdf_models$future

  all_bins <- 1:nbins
  selected_indices_all <- lapply(seq_along(lon), function(i) lapply(seq_along(lat), function(j) all_bins))

  h_dist_present <- compute_partial_hdist(pdf_ref_present, pdf_models_present, selected_indices_all)
  h_dist_future  <- compute_partial_hdist(pdf_ref_future, pdf_models_future, selected_indices_all)

  hist(h_dist_present, main = "Complete Hellinger Distance (Present)")
  hist(h_dist_future, main = "Complete Hellinger Distance (Future)")

  MMM_present <- apply(pdf_models_present, c(1, 2, 3), mean)
  MMM_future <- apply(pdf_models_future, c(1, 2, 3), mean)

  MMM_hdist_present <- matrix(NA, length(lon), length(lat))
  MMM_hdist_future <- matrix(NA, length(lon), length(lat))
  for (i in seq_along(lon)) {
    for (j in seq_along(lat)) {
      MMM_hdist_present[i, j] <- sqrt(sum((sqrt(MMM_present[i, j, ]) - sqrt(pdf_ref_present[i, j, ]))^2)) / sqrt(2)
      MMM_hdist_future[i, j] <- sqrt(sum((sqrt(MMM_future[i, j, ]) - sqrt(pdf_ref_future[i, j, ]))^2)) / sqrt(2)
    }
  }

  ldr_indices <- lapply(seq_len(lon_size), function(i) {
    lapply(seq_len(lat_size), function(j) {
      central_indices <- select_hdr_indices(pdf_ref_future[i, j, ], tau = 0.10)
      setdiff(seq_len(nbins), central_indices)
    })
  })

  partial_hdist_future <- compute_partial_hdist(pdf_ref_future, pdf_models_future, ldr_indices)

  for (smooth_cost in smooth_costs) {
    result <- tryCatch({
      GraphCutHellinger_nD_lat(
        pdf_models_future = pdf_models_future,
        h_dist = h_dist_present,
        weight_data = 1,
        weight_smooth = smooth_cost,
        nBins = nbins,
        lat = lat,
        seed = 1,
        verbose = TRUE,
        rebuild = TRUE
      )
    }, error = function(e) {
      cat("Error with smooth_cost =", smooth_cost, ":", e$message, "\n")
      NULL
    })

    if (!is.null(result)) {
      label_name <- as.character(smooth_cost)

      hdist_present <- matrix(NA, length(lon), length(lat))
      hdist_future <- matrix(NA, length(lon), length(lat))
      partial_h <- matrix(NA, length(lon), length(lat))

      for (l in seq_along(model_names)) {
        islabel <- which(result$label_attribution == l)
        hdist_present[islabel] <- h_dist_present[,,l][islabel]
        hdist_future[islabel] <- h_dist_future[,,l][islabel]
        partial_h[islabel] <- partial_hdist_future[,,l][islabel]
      }

      GC_result_list[[label_name]][[reference_name]] <- result
      GC_hdist_present_list[[label_name]][[reference_name]] <- hdist_present
      GC_hdist_future_list[[label_name]][[reference_name]] <- hdist_future
      GC_partial_hdist_future_list[[label_name]][[reference_name]] <- partial_h
    }
    gc()
  }

  MMM_partial_hdist_future <- matrix(NA, lon_size, lat_size)
  for (i in seq_len(lon_size)) {
    for (j in seq_len(lat_size)) {
      selected_bins <- ldr_indices[[i]][[j]]
      if (length(selected_bins) > 0) {
        MMM_partial_hdist_future[i, j] <- sqrt(
          sum((sqrt(MMM_future[i, j, selected_bins]) - sqrt(pdf_ref_future[i, j, selected_bins]))^2)
        ) / sqrt(2)
      }
    }
  }

  cat("Mean Partial Hellinger Distance for MMM (future):", mean(MMM_partial_hdist_future, na.rm = TRUE), "\n")

  MMM_hdist_present_list[[reference_name]] <- MMM_hdist_present
  MMM_hdist_future_list[[reference_name]] <- MMM_hdist_future
  MMM_partial_hdist_future_list[[reference_name]] <- MMM_partial_hdist_future
  MMM_present_list[[reference_name]] <- MMM_present
  MMM_future_list[[reference_name]] <- MMM_future


  # Sauvegarde de sécurité après chaque modèle traité
  ref_short <- gsub("[^A-Za-z0-9]", "", reference_name)
  checkpoint_file <- paste0("checkpoints/PME_", ref_short, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")

  # Crée le dossier 'checkpoints' s'il n'existe pas
  if (!dir.exists("checkpoints")) dir.create("checkpoints")

  saveRDS(
    list(
      GC_result_list = GC_result_list,
      GC_hdist_future_list = GC_hdist_future_list,
      GC_hdist_present_list = GC_hdist_present_list,
      GC_partial_hdist_future_list = GC_partial_hdist_future_list,
      MMM_hdist_future_list = MMM_hdist_future_list,
      MMM_hdist_present_list = MMM_hdist_present_list,
      MMM_partial_hdist_future_list = MMM_partial_hdist_future_list,
      MMM_present_list = MMM_present_list,
      MMM_future_list = MMM_future_list
    ),
    file = checkpoint_file,
    compress = FALSE
  )
  cat("Checkpoint saved to", checkpoint_file, "\n")

}



# Get the current date and time
current_time <- Sys.time()
# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")
# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "_my_workspace_PME_bias_corrected_22models_full.RData")
# Save the workspace using the generated filename
save.image(file = filename, compress = FALSE)



# Boxplot of the H dist projection by ref
# Get the name of the first reference model
ref_name <- names(GC06_hdist_future_list)[1]

# Extract the corresponding Hellinger matrices
gc_h <- GC06_hdist_future_list[[ref_name]]
mmm_h <- MMM_hdist_future_list[[ref_name]]

# Flatten the matrices into vectors
gc_values <- as.vector(gc_h)
mmm_values <- as.vector(mmm_h)

# Combine into a data frame for ggplot
df <- data.frame(
  Hellinger = c(gc_values, mmm_values),
  Method = factor(rep(c("GraphCut", "MMM"), each = length(gc_values)))
)

# Remove NAs if any
df <- na.omit(df)

# Plot
p <- ggplot(df, aes(x = Method, y = Hellinger, fill = Method)) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.3) +
  scale_fill_manual(values = c("GraphCut" = "#0072B2", "MMM" = "#D55E00")) +
  labs(
    title = paste("Hellinger Distance (Future) - Reference:", ref_name),  # <-- fixed here
    y = "Hellinger Distance",
    x = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

print(p)


# Violin plot of the H dist by ref
{
  # Boxplot of the H dist projection by ref
  # Get the name of the first reference model
  ref_name <- names(GC06_hdist_future_list)[1]

  # Extract the corresponding Hellinger matrices
  gc_h <- GC06_hdist_future_list[[ref_name]]
  mmm_h <- MMM_hdist_future_list[[ref_name]]

  # Flatten the matrices into vectors
  gc_values <- as.vector(gc_h)
  mmm_values <- as.vector(mmm_h)

  # Combine into a data frame for ggplot

  # Rebuild df_all with full Hellinger distances for all references
  df_all <- do.call(rbind, lapply(names(GC06_hdist_future_list), function(ref_name) {
    gc_vals <- as.vector(GC06_hdist_future_list[[ref_name]])
    mmm_vals <- as.vector(MMM_hdist_future_list[[ref_name]])

    data.frame(
      Hellinger = c(gc_vals, mmm_vals),
      Method = factor(rep(c("GraphCut", "MMM"), each = length(gc_vals))),
      Reference = ref_name
    )
  }))
  df_all <- na.omit(df_all)
  library(ggplot2)
  library(dplyr)

  df_all$Reference <- factor(df_all$Reference, levels = unique(df_all$Reference))

  ggplot(df_all, aes(x = Reference, y = Hellinger, fill = Method)) +
    geom_violin(position = position_dodge(width = 0.6),
                width = 2,      # Plus large horizontalement
                adjust = 1,     # Lissage plus doux (tu peux tester 1, 1.5, 2)
                alpha = 0.79) +
    scale_fill_manual(values = c("GraphCut" = "#0072B2", "MMM" = "#D55E00")) +
    labs(
      title = "Distribution of Hellinger Distance (Future) by Reference Model",
      x = "Reference Model",
      y = "Hellinger Distance"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
      legend.position = "top",
      plot.title = element_text(size = 16, face = "bold")
    )

}

# Violin plot of the partial H dist by ref
{

  # Boxplot of the H dist projection by ref
  # Get the name of the first reference model
  ref_name <- names(GC06_hdist_future_list)[1]

  # Extract the corresponding Hellinger matrices
  gc_h <- GC06_partial_hdist_future_list[[ref_name]]
  mmm_h <- MMM_partial_hdist_future_list[[ref_name]]

  # Flatten the matrices into vectors
  gc_values <- as.vector(gc_h)
  mmm_values <- as.vector(mmm_h)

  # Rebuild df_all with partial Hellinger distances
  df_all <- do.call(rbind, lapply(names(GC06_partial_hdist_future_list), function(ref_name) {
    gc_vals <- as.vector(GC06_partial_hdist_future_list[[ref_name]])
    mmm_vals <- as.vector(MMM_partial_hdist_future_list[[ref_name]])

    data.frame(
      Hellinger = c(gc_vals, mmm_vals),
      Method = factor(rep(c("GraphCut", "MMM"), each = length(gc_vals))),
      Reference = ref_name
    )
  }))
  df_all <- na.omit(df_all)

  library(ggplot2)
  library(dplyr)

  df_all$Reference <- factor(df_all$Reference, levels = unique(df_all$Reference))

  ggplot(df_all, aes(x = Reference, y = Hellinger, fill = Method)) +
    geom_violin(position = position_dodge(width = 0.6),
                width = 2,      # Plus large horizontalement
                adjust = 1,     # Lissage plus doux (tu peux tester 1, 1.5, 2)
                alpha = 0.79) +
    scale_fill_manual(values = c("GraphCut" = "#0072B2", "MMM" = "#D55E00")) +
    labs(
      title = "Distribution of Hellinger Distance (Future) by Reference Model",
      x = "Reference Model",
      y = "Hellinger Distance"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
      legend.position = "top",
      plot.title = element_text(size = 16, face = "bold")
    )

}

# Aggregated violin by method
{
  # Rebuild df_all with full Hellinger distances for all references
  df_all <- do.call(rbind, lapply(names(GC06_hdist_future_list), function(ref_name) {
    gc_vals <- as.vector(GC06_hdist_future_list[[ref_name]])
    mmm_vals <- as.vector(MMM_hdist_future_list[[ref_name]])

    data.frame(
      Hellinger = c(gc_vals, mmm_vals),
      Method = factor(rep(c("GraphCut", "MMM"), each = length(gc_vals))),
      Reference = ref_name
    )
  }))
  df_all <- na.omit(df_all)
  library(ggplot2)
  library(dplyr)

  # Flatten both lists into a long data.frame with method and value
  gc_data <- do.call(c, lapply(GC06_hdist_future_list, as.vector))
  mmm_data <- do.call(c, lapply(MMM_hdist_future_list, as.vector))

  # Remove NA values
  gc_data <- gc_data[!is.na(gc_data)]
  mmm_data <- mmm_data[!is.na(mmm_data)]

  # Combine into a data.frame
  df_all_methods <- data.frame(
    Method = factor(rep(c("GraphCut", "MMM"), times = c(length(gc_data), length(mmm_data))),
                    levels = c("GraphCut", "MMM")),
    Hellinger = c(gc_data, mmm_data)
  )

  # Plot
  ggplot(df_all_methods, aes(x = Method, y = Hellinger, fill = Method)) +
    geom_violin(scale = "area", trim = TRUE, adjust = 1.5) +
    scale_fill_manual(values = c("GraphCut" = "#1f78b4", "MMM" = "#e66101")) +
    theme_minimal(base_size = 14) +
    labs(
      title = "Distribution of Hellinger Distance (Future)",
      subtitle = "Aggregated across all reference models",
      y = "Hellinger Distance",
      x = NULL
    )

}


# Violin plot of the H dist by ref multiple smooth cost
{
  library(ggplot2)
  library(dplyr)

  # Crée le dossier s'il n'existe pas
  if (!dir.exists("figure")) dir.create("figure")

  # Liste des références disponibles
  references_done <- names(MMM_hdist_future_list)

  # Boucle sur chaque référence
  for (ref_name in references_done) {

    # Créer un dataframe pour cette référence
    mmm_vals <- as.vector(MMM_hdist_future_list[[ref_name]])
    df_ref <- data.frame(
      Hellinger = mmm_vals,
      Method = "MMM",
      SmoothCost = "MMM",
      Reference = ref_name
    )

    for (cost in c("0.05", "0.1", "0.6", "1", "2")) {
      if (!is.null(GC_hdist_future_list[[cost]][[ref_name]])) {
        gc_vals <- as.vector(GC_hdist_future_list[[cost]][[ref_name]])
        df_ref <- rbind(df_ref, data.frame(
          Hellinger = gc_vals,
          Method = "GraphCut",
          SmoothCost = cost,
          Reference = ref_name
        ))
      }
    }

    df_ref <- na.omit(df_ref)
    df_ref$MethodLabel <- ifelse(df_ref$Method == "MMM", "MMM", paste0("GC-", df_ref$SmoothCost))
    df_ref$MethodLabel <- factor(df_ref$MethodLabel, levels = c("GC-0.05", "GC-0.1", "GC-0.6", "GC-1", "GC-2", "MMM"))

    # Générer le plot
    p_ref <- ggplot(df_ref, aes(x = MethodLabel, y = Hellinger, fill = MethodLabel)) +
      geom_violin(scale = "area", adjust = 1.2, width = 0.7, alpha = 0.85) +
      # Ajout moyenne (point noir)
      stat_summary(fun = mean, geom = "point", shape = 20, size = 2.2, color = "black", position = position_dodge(width = 0.7))+
      # Ajout médiane (barre rouge)
      stat_summary(fun = median, geom = "crossbar", width = 0.4, color = "red", fatten = 1, position = position_dodge(width = 0.7))+
      scale_fill_manual(values = c(
        "GC-0.05" = "#1b9e77",
        "GC-0.1"  = "#d95f02",
        "GC-0.6"  = "#7570b3",
        "GC-1"    = "#e7298a",
        "GC-2"    = "#66a61e",
        "MMM"     = "#e6ab02"
      )) +
      labs(
        title = paste("Hellinger Distance (Future) - Reference:", ref_name),
        x = "Method",
        y = "Hellinger Distance",
        fill = "Method"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(size = 10),
        plot.title = element_text(size = 16, face = "bold"),
        legend.position = "none"
      )

    # Nettoyer le nom du fichier
    clean_name <- gsub("[^A-Za-z0-9]", "", ref_name)
    file_base <- paste0("figure/Violin_Hdist_BC_22model_", clean_name)

    # Sauvegarde en PDF et PNG
    ggsave(paste0(file_base, ".pdf"), plot = p_ref, width = 20, height = 15, units = "cm", dpi = 300)
    ggsave(paste0(file_base, ".png"), plot = p_ref, width = 20, height = 15, units = "cm", dpi = 300)

    cat("Saved plot for reference:", ref_name, "\n")
  }

}


#Violin partial H by ref
{
  library(ggplot2)
  library(dplyr)

  # Crée le dossier 'figure' s'il n'existe pas
  if (!dir.exists("figure")) dir.create("figure")

  # Liste des références disponibles
  references_done <- names(MMM_partial_hdist_future_list)

  for (ref_name in references_done) {

    # MMM
    mmm_vals <- as.vector(MMM_partial_hdist_future_list[[ref_name]])
    df_ref <- data.frame(
      Hellinger = mmm_vals,
      Method = "MMM",
      SmoothCost = "MMM",
      Reference = ref_name
    )

    # Tous les GraphCut
    for (cost in c("0.05", "0.1", "0.6", "1", "2")) {
      if (!is.null(GC_partial_hdist_future_list[[cost]][[ref_name]])) {
        gc_vals <- as.vector(GC_partial_hdist_future_list[[cost]][[ref_name]])
        df_ref <- rbind(df_ref, data.frame(
          Hellinger = gc_vals,
          Method = "GraphCut",
          SmoothCost = cost,
          Reference = ref_name
        ))
      }
    }

    df_ref <- na.omit(df_ref)
    df_ref$MethodLabel <- ifelse(df_ref$Method == "MMM", "MMM", paste0("GC-", df_ref$SmoothCost))
    df_ref$MethodLabel <- factor(df_ref$MethodLabel, levels = c("GC-0.05", "GC-0.1", "GC-0.6", "GC-1", "GC-2", "MMM"))

    # Plot
    p_ref <- ggplot(df_ref, aes(x = MethodLabel, y = Hellinger, fill = MethodLabel)) +
      geom_violin(scale = "area", adjust = 1.2, width = 0.7, alpha = 0.85) +
      scale_fill_manual(values = c(
        "GC-0.05" = "#1b9e77",
        "GC-0.1"  = "#d95f02",
        "GC-0.6"  = "#7570b3",
        "GC-1"    = "#e7298a",
        "GC-2"    = "#66a61e",
        "MMM"     = "#e6ab02"
      )) +
      labs(
        title = paste("Partial Hellinger Distance (Future) - Reference:", ref_name),
        x = "Method",
        y = "Partial Hellinger Distance",
        fill = "Method"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(size = 10),
        plot.title = element_text(size = 16, face = "bold"),
        legend.position = "none"
      )

    # Nettoyer le nom du fichier
    clean_name <- gsub("[^A-Za-z0-9]", "", ref_name)
    file_base <- paste0("figure/Violin_Hdist_BC_22model_", clean_name)

    # Sauvegarde
    ggsave(paste0(file_base, "_partial.pdf"), plot = p_ref, width = 20, height = 15, units = "cm", dpi = 300)
    ggsave(paste0(file_base, "_partial.png"), plot = p_ref, width = 20, height = 15, units = "cm", dpi = 300)

    cat("Saved partial Hellinger violin plot for reference:", ref_name, "\n")
  }

}

# Violins aggragated by method
{
  library(ggplot2)
  library(dplyr)

  # Créer le dossier si nécessaire
  if (!dir.exists("figure")) dir.create("figure")

  # Construire un data frame agrégé : toutes méthodes, toutes références
  df_all <- do.call(rbind, lapply(names(MMM_hdist_future_list), function(ref_name) {
    mmm_vals <- as.vector(MMM_hdist_future_list[[ref_name]])
    df <- data.frame(
      Hellinger = mmm_vals,
      Method = "MMM",
      SmoothCost = "MMM",
      Reference = ref_name
    )

    for (cost in c("0.05", "0.1", "0.6", "1", "2")) {
      if (!is.null(GC_hdist_future_list[[cost]][[ref_name]])) {
        gc_vals <- as.vector(GC_hdist_future_list[[cost]][[ref_name]])
        df <- rbind(df, data.frame(
          Hellinger = gc_vals,
          Method = "GraphCut",
          SmoothCost = cost,
          Reference = ref_name
        ))
      }
    }

    return(df)
  }))

  df_all <- na.omit(df_all)
  df_all$MethodLabel <- ifelse(df_all$Method == "MMM", "MMM", paste0("GC-", df_all$SmoothCost))
  df_all$MethodLabel <- factor(df_all$MethodLabel, levels = c("GC-0.05", "GC-0.1", "GC-0.6", "GC-1", "GC-2", "MMM"))

  # Créer le plot agrégé
  p_agg <- ggplot(df_all, aes(x = MethodLabel, y = Hellinger, fill = MethodLabel)) +
    geom_violin(scale = "area", trim = TRUE, adjust = 1.5, alpha = 0.85, width = 0.7) +

    # Ajout moyenne (point bleu foncé)
    stat_summary(fun = mean, geom = "point", shape = 20, size = 2.2, color = "black", position = position_dodge(width = 0.7)) +

    # Ajout médiane (barre rouge)
    stat_summary(fun = median, geom = "crossbar", width = 0.4, color = "red", fatten = 1, position = position_dodge(width = 0.7)) +

    scale_fill_manual(values = c(
      "GC-0.05" = "#1b9e77",
      "GC-0.1"  = "#d95f02",
      "GC-0.6"  = "#7570b3",
      "GC-1"    = "#e7298a",
      "GC-2"    = "#66a61e",
      "MMM"     = "#e6ab02"
    )) +
    labs(
      title = "Distribution of Hellinger Distance (Future)",
      subtitle = "Aggregated across all reference models",
      x = "Method",
      y = "Hellinger Distance",
      fill = "Method"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 11),
      plot.title = element_text(size = 16, face = "bold"),
      legend.position = "none"
    )

  # Sauvegarde
  file_base <- "figure/Aggregated_Hdist_BC_22models"
  ggsave(paste0(file_base, ".pdf"), plot = p_agg, width = 20, height = 15, units = "cm", dpi = 300)
  ggsave(paste0(file_base, ".png"), plot = p_agg, width = 20, height = 15, units = "cm", dpi = 300)

  cat("Saved aggregated Hellinger distance plot to", file_base, "\n")

}


# Aggragated partial H
{

  library(ggplot2)
  library(dplyr)

  # Créer le dossier si nécessaire
  if (!dir.exists("figure")) dir.create("figure")

  # Construire un data frame agrégé : toutes méthodes, toutes références
  df_all <- do.call(rbind, lapply(names(MMM_partial_hdist_future_list), function(ref_name) {
    mmm_vals <- as.vector(MMM_partial_hdist_future_list[[ref_name]])
    df <- data.frame(
      Hellinger = mmm_vals,
      Method = "MMM",
      SmoothCost = "MMM",
      Reference = ref_name
    )

    for (cost in c("0.05", "0.1", "0.6", "1", "2")) {
      if (!is.null(GC_partial_hdist_future_list[[cost]][[ref_name]])) {
        gc_vals <- as.vector(GC_partial_hdist_future_list[[cost]][[ref_name]])
        df <- rbind(df, data.frame(
          Hellinger = gc_vals,
          Method = "GraphCut",
          SmoothCost = cost,
          Reference = ref_name
        ))
      }
    }

    return(df)
  }))

  df_all <- na.omit(df_all)
  df_all$MethodLabel <- ifelse(df_all$Method == "MMM", "MMM", paste0("GC-", df_all$SmoothCost))
  df_all$MethodLabel <- factor(df_all$MethodLabel, levels = c("GC-0.05", "GC-0.1", "GC-0.6", "GC-1", "GC-2", "MMM"))

  # Créer le plot agrégé
  p_agg <- ggplot(df_all, aes(x = MethodLabel, y = Hellinger, fill = MethodLabel)) +
    geom_violin(scale = "area", trim = TRUE, adjust = 1.5, alpha = 0.85, width = 0.7) +

    # Ajout moyenne (point noir)
    stat_summary(fun = mean, geom = "point", shape = 20, size = 2.2, color = "black", position = position_dodge(width = 0.7)) +

    # Ajout médiane (barre rouge)
    stat_summary(fun = median, geom = "crossbar", width = 0.4, color = "red", fatten = 1, position = position_dodge(width = 0.7)) +

    scale_fill_manual(values = c(
      "GC-0.05" = "#1b9e77",
      "GC-0.1"  = "#d95f02",
      "GC-0.6"  = "#7570b3",
      "GC-1"    = "#e7298a",
      "GC-2"    = "#66a61e",
      "MMM"     = "#e6ab02"
    )) +
    labs(
      title = "Partial Hellinger Distance (Future)",
      subtitle = "Aggregated across all reference models",
      x = "Method",
      y = "Partial Hellinger Distance",
      fill = "Method"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 11),
      plot.title = element_text(size = 16, face = "bold"),
      legend.position = "none"
    )

  # Sauvegarde
  file_base <- "figure/Aggregated_HdistPartial_BC_22models"
  ggsave(paste0(file_base, ".pdf"), plot = p_agg, width = 20, height = 15, units = "cm", dpi = 300)
  ggsave(paste0(file_base, ".png"), plot = p_agg, width = 20, height = 15, units = "cm", dpi = 300)

  cat("Saved aggregated partial Hellinger distance plot to", file_base, "\n")


}

# violin panel with 22 ref : H dist
{
  library(ggplot2)
  library(dplyr)
  library(patchwork)

  # Liste des références disponibles
  references_done <- names(MMM_hdist_future_list)

  # Créer tous les plots dans une liste
  plot_list <- list()

  for (ref_name in references_done) {

    mmm_vals <- as.vector(MMM_hdist_future_list[[ref_name]])
    df_ref <- data.frame(
      Hellinger = mmm_vals,
      Method = "MMM",
      SmoothCost = "MMM",
      Reference = ref_name
    )

    for (cost in c("0.05", "0.1", "0.6", "1", "2")) {
      if (!is.null(GC_hdist_future_list[[cost]][[ref_name]])) {
        gc_vals <- as.vector(GC_hdist_future_list[[cost]][[ref_name]])
        df_ref <- rbind(df_ref, data.frame(
          Hellinger = gc_vals,
          Method = "GraphCut",
          SmoothCost = cost,
          Reference = ref_name
        ))
      }
    }

    df_ref <- na.omit(df_ref)
    df_ref$MethodLabel <- ifelse(df_ref$Method == "MMM", "MMM", paste0("GC-", df_ref$SmoothCost))
    df_ref$MethodLabel <- factor(df_ref$MethodLabel, levels = c("GC-0.05", "GC-0.1", "GC-0.6", "GC-1", "GC-2", "MMM"))

    # Create the individual plot
    p <- ggplot(df_ref, aes(x = MethodLabel, y = Hellinger, fill = MethodLabel)) +
      geom_violin(scale = "area", adjust = 1.2, width = 0.7, alpha = 0.85) +
      stat_summary(fun = mean, geom = "point", shape = 20, size = 2.2, color = "black", position = position_dodge(width = 0.7)) +
      stat_summary(fun = median, geom = "crossbar", width = 0.4, color = "red", fatten = 1, position = position_dodge(width = 0.7)) +
      scale_fill_manual(values = c(
        "GC-0.05" = "#1b9e77",
        "GC-0.1"  = "#d95f02",
        "GC-0.6"  = "#7570b3",
        "GC-1"    = "#e7298a",
        "GC-2"    = "#66a61e",
        "MMM"     = "#e6ab02"
      )) +
      labs(
        title = ref_name,
        x = NULL,
        y = NULL
      ) +
      theme_minimal(base_size = 10) +
      theme(
        axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(size = 10, face = "bold"),
        legend.position = "none"
      )

    plot_list[[ref_name]] <- p
  }

  # Combine all plots with patchwork
  big_plot <- wrap_plots(plot_list, ncol = 8) +
    plot_annotation(
      title = "Hellinger Distance (Future) by Reference and Method",
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
      )
    )

  # Save in high-resolution (4K scale, roughly)
  ggsave("figure/AllReferences_Hellinger_Grid_22models.png", plot = big_plot,
         width = 3840/96, height = 2160/96, dpi = 300, units = "in")  # 4K: 3840×2160 pixels
  ggsave("figure/AllReferences_Hellinger_Grid_22models.pdf", plot = big_plot,
         width = 3840/96, height = 2160/96, dpi = 300, units = "in")

  cat("✅ Multi-panel violin plot saved as 4K image and PDF.\n")

}

{
  library(ggplot2)
  library(dplyr)
  library(patchwork)

  # Liste des références disponibles
  references_done <- names(MMM_partial_hdist_future_list)

  # Créer tous les plots dans une liste
  plot_list <- list()

  for (ref_name in references_done) {

    mmm_vals <- as.vector(MMM_partial_hdist_future_list[[ref_name]])
    df_ref <- data.frame(
      Hellinger = mmm_vals,
      Method = "MMM",
      SmoothCost = "MMM",
      Reference = ref_name
    )

    for (cost in c("0.05", "0.1", "0.6", "1", "2")) {
      if (!is.null(GC_partial_hdist_future_list[[cost]][[ref_name]])) {
        gc_vals <- as.vector(GC_partial_hdist_future_list[[cost]][[ref_name]])
        df_ref <- rbind(df_ref, data.frame(
          Hellinger = gc_vals,
          Method = "GraphCut",
          SmoothCost = cost,
          Reference = ref_name
        ))
      }
    }

    df_ref <- na.omit(df_ref)
    df_ref$MethodLabel <- ifelse(df_ref$Method == "MMM", "MMM", paste0("GC-", df_ref$SmoothCost))
    df_ref$MethodLabel <- factor(df_ref$MethodLabel, levels = c("GC-0.05", "GC-0.1", "GC-0.6", "GC-1", "GC-2", "MMM"))

    # Créer le plot individuel
    p <- ggplot(df_ref, aes(x = MethodLabel, y = Hellinger, fill = MethodLabel)) +
      geom_violin(scale = "area", adjust = 1.2, width = 0.7, alpha = 0.85) +
      stat_summary(fun = mean, geom = "point", shape = 20, size = 2.2, color = "black", position = position_dodge(width = 0.7)) +
      stat_summary(fun = median, geom = "crossbar", width = 0.4, color = "red", fatten = 1, position = position_dodge(width = 0.7)) +
      scale_fill_manual(values = c(
        "GC-0.05" = "#1b9e77",
        "GC-0.1"  = "#d95f02",
        "GC-0.6"  = "#7570b3",
        "GC-1"    = "#e7298a",
        "GC-2"    = "#66a61e",
        "MMM"     = "#e6ab02"
      )) +
      labs(
        title = ref_name,
        x = NULL,
        y = NULL
      ) +
      theme_minimal(base_size = 10) +
      theme(
        axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(size = 10, face = "bold"),
        legend.position = "none"
      )

    plot_list[[ref_name]] <- p
  }

  # Combinaison avec patchwork
  big_plot <- wrap_plots(plot_list, ncol = 8) +
    plot_annotation(
      title = "Partial Hellinger Distance (Future) by Reference and Method",
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
      )
    )

  # Sauvegarde en 4K
  file_base <- "figure/AllReferences_HellingerPartial_Grid_22models"
  ggsave(paste0(file_base, ".png"), plot = big_plot,
         width = 3840/96, height = 2160/96, dpi = 300, units = "in")
  ggsave(paste0(file_base, ".pdf"), plot = big_plot,
         width = 3840/96, height = 2160/96, dpi = 300, units = "in")

  cat("✅ Multi-panel partial Hellinger plot saved as 4K image and PDF.\n")

}

# Summary table
{
  # Résumé des statistiques Hellinger par méthode et par référence
  summary_stats <- df_all %>%
    group_by(Reference, MethodLabel) %>%
    summarise(
      Mean = mean(Hellinger, na.rm = TRUE),
      Median = median(Hellinger, na.rm = TRUE),
      .groups = "drop"
    )

  # Ajouter la ligne globale (toutes références confondues)
  global_summary <- df_all %>%
    group_by(MethodLabel) %>%
    summarise(
      Mean = mean(Hellinger, na.rm = TRUE),
      Median = median(Hellinger, na.rm = TRUE)
    ) %>%
    mutate(Reference = "ALL")

  # Combiner
  summary_stats <- bind_rows(summary_stats, global_summary)

  # Afficher le tableau
  print(summary_stats)

  # Optionnel : export CSV
  # write.csv(summary_stats, "summary_Hellinger_stats.csv", row.names = FALSE)

}

{
  library(ggplot2)
  library(dplyr)
  library(tidyr)

  # Reordonner les références avec "ALL" à la fin
  summary_stats <- summary_stats %>%
    mutate(
      Reference = factor(Reference, levels = c(setdiff(unique(Reference), "ALL"), "ALL")),
      MethodLabel = factor(MethodLabel, levels = c("GC-0.05", "GC-0.1", "GC-0.6", "GC-1", "GC-2", "MMM"))
    )

  # Identifier la meilleure méthode par référence (plus petite moyenne)
  summary_stats_highlight <- summary_stats %>%
    group_by(Reference) %>%
    mutate(
      IsBest = Mean == min(Mean, na.rm = TRUE),
      FontFace = ifelse(IsBest, "bold", "plain")
    ) %>%
    ungroup()

  # Heatmap avec mise en évidence
  p_mean <- ggplot(summary_stats_highlight, aes(x = MethodLabel, y = Reference, fill = Mean)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.3f", Mean), fontface = FontFace), size = 3, color = "black") +
    scale_fill_viridis_c(option = "C", name = "Mean\nHellinger") +
    labs(
      title = "Mean Hellinger Distance by Method and Reference",
      subtitle = "Best per row is shown in bold",
      x = "Method",
      y = "Reference"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )

  # Affichage
  print(p_mean)

  # Sauvegarde
  ggsave("figure/Summary_Hellinger_Heatmap_Mean_Highlighted.pdf", plot = p_mean, width = 20, height = 15, units = "cm", dpi = 300)
  ggsave("figure/Summary_Hellinger_Heatmap_Mean_Highlighted.png", plot = p_mean, width = 20, height = 15, units = "cm", dpi = 300)
}


{
  library(ggplot2)
  library(dplyr)
  library(tidyr)

  # Identifier la méthode gagnante par référence (plus petite médiane)
  summary_stats_highlight <- summary_stats %>%
    group_by(Reference) %>%
    mutate(
      IsBest = Median == min(Median, na.rm = TRUE),
      FontFace = ifelse(IsBest, "bold", "plain")
    ) %>%
    ungroup()

  # Heatmap avec texte noir et valeur gagnante en gras
  p_median <- ggplot(summary_stats_highlight, aes(x = MethodLabel, y = Reference, fill = Median)) +
    geom_tile(color = "white") +
    geom_text(
      aes(label = sprintf("%.3f", Median), fontface = FontFace),
      size = 3, color = "black"
    ) +
    scale_fill_viridis_c(option = "C", name = "Median\nHellinger") +
    labs(
      title = "Median Hellinger Distance by Method and Reference",
      subtitle = "Best per row is shown in bold",
      x = "Method",
      y = "Reference"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )

  # Affichage + sauvegarde
  print(p_median)
  ggsave("figure/Summary_Hellinger_Heatmap_Median_Highlighted.pdf", plot = p_median,
         width = 20, height = 15, units = "cm", dpi = 300)
  ggsave("figure/Summary_Hellinger_Heatmap_Median_Highlighted.png", plot = p_median,
         width = 20, height = 15, units = "cm", dpi = 300)
}

