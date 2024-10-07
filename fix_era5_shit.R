tmp <- OpenAndHist2D_range('ERA5', variables, year_present, range_var_final)

kde_ref <- tmp[[1]]


tmp <- OpenAndHist2D_range('ERA5', variables, year_future , range_var_final)

kde_ref_future <- tmp[[1]]

rm(pdf_matrix)
rm(tmp)


tmp <- OpenAndAverageCMIP6(
  reference_name, variables, year_present, year_future
)
reference_list <- tmp[[1]]
reference_matrix <- tmp[[2]]
rm(tmp)

models_matrix_nrm <- list()
models_matrix_nrm <- NormalizeVariables(models_matrix, variables, 'StdSc')

reference_matrix_nrm <- list()
reference_matrix_nrm <- NormalizeVariables(reference_matrix, variables, 'StdSc')


weight_data <- 1
weight_smooth <- 1


GC_result_hellinger_test_new_era5_5 <- list()
GC_result_hellinger_test_new_era5_5 <- tryCatch({
  # Call the GraphCutHellinger2D_new2 function and store the result
  GraphCutHellinger2D_new2(pdf_ref = kde_ref,
                           kde_models = kde_models,
                           pdf_models_future = kde_models_future,
                           h_dist = h_dist,
                           weight_data = weight_data,
                           weight_smooth = weight_smooth,
                           nBins = nbins1d^2,
                           verbose = TRUE,
                           rebuild = FALSE)
}, error = function(e) {
  # If an error occurs, save the error message and return NULL for this iteration
  error_list[[i]] <- paste("Error in iteration", i, ":", e$message)
  NULL  # Returning NULL to indicate failure for this iteration
})



GC_hellinger_projections <- list()
j <- 1
for(var in variables){
  for(l in 0:(length(model_names))){
    islabel <- which(GC_result_hellinger_test_new_era5_5$label_attribution == l)
    GC_hellinger_projections[[var]][islabel] <- models_matrix$future[,,(l),j][islabel]
  }
  j <- j + 1
}
GC_hellinger_projections$tas <- matrix(GC_hellinger_projections$tas, nrow = 360)
GC_hellinger_projections$pr <- matrix(GC_hellinger_projections$pr * 86400, nrow = 360)





h_dist_map <- array(NA, dim = dim(GC_result_hellinger_test_new_era5_5$label_attribution))

for(j in seq_along(variables)){
  for(l in 0:(length(model_names))){
    islabel <- which(GC_result_hellinger_test_new_era5_5$label_attribution == l)
    h_dist_map[islabel] <- h_dist_future[,,(l)][islabel]
  }
}

test_df <- melt(h_dist_map, c("lon", "lat"), value.name = "Bias")



p5 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias))+
  labs(subtitle = 'Projection period : 1999 - 2022')+
  ggtitle(paste0('GraphCut Hellinger', ': Mean Hellinger distance = ', round(mean(h_dist_map), 2)))+
  scale_fill_gradient(low = "white", high = "#015a8c", limits = c(0.1, 0.70), oob = scales::squish)+
  borders("world2", colour = 'black', lwd = 0.12) +
  scale_x_continuous(, expand = c(0, 0)) +
  scale_y_continuous(, expand = c(0,0))+
  theme(legend.position = 'bottom')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())+
  xlab('Longitude')+
  ylab('Latitude') +
  labs(fill='Hellinger \nDistance')+
  theme_bw()+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=16), #change legend title font sizen
        legend.text = element_text(size=12))+ #change legend text font size
  theme(plot.title = element_text(size=24),
        plot.subtitle = element_text(size = 20,hjust=0.5),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),)+
  easy_center_title()
p5

name <- paste0('figure/H_dist_future_GC_hellinger_new_26models5')
ggsave(paste0(name, '.pdf'), plot = p5, width = 35, height = 25, units = "cm", dpi = 300)
ggsave(paste0(name, '.png'), plot = p5, width = 35, height = 25, units = "cm", dpi = 300)



# Map of the bias that works
# Creating the dataframe
# For tas
bias_tmp <- GC_hellinger_projections$tas - reference_list$future$tas[[1]]

test_df <- melt(bias_tmp, c("lon", "lat"), value.name = "Bias")

# Settings the limits
# limit <- max(abs(c(min(test_df$Bias), max(test_df$Bias))))
# limit <- round_any(limit, 10, floor)

limit <- 8
limits <- c(-8, 8)
v_limits <- as.numeric(format(seq(-limit, limit, len=6), digits = 3))

test_df$Bias[test_df$Bias < -limit] <- -limit
test_df$Bias[test_df$Bias > limit] <- limit

p <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias),)+
  ggtitle(paste0('GraphCut Hellinger'))+
  labs(subtitle = paste0('Average Bias : ', round(mean(abs(bias_tmp)), 2)))+
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),
                       breaks = v_limits,
                       limits = limits)+
  # guide = guide_colourbar(
  #   barwidth = 0.4,
  #   ticks.colour = 'black',
  #   ticks.linewidth = 1,
  #   frame.colour = 'black',
  #   draw.ulim = TRUE,
  #   draw.llim = TRUE),)
  borders("world2", colour = 'black', lwd = 0.12) +
  scale_x_continuous(, expand = c(0, 0)) +
  scale_y_continuous(, expand = c(0,0))+
  theme(legend.position = 'bottom')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())+
  xlab('Longitude')+
  ylab('Latitude') +
  labs(fill='Bias [K]')+
  theme_bw()+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=12))+ #change legend text font size
  theme(plot.title = element_text(size=24),
        plot.subtitle = element_text(size = 20,hjust=0.5),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),)+
  easy_center_title()
p
mean(abs(bias_tmp))

name <- paste0('figure/GC_Hellinger_bias_tas_26models_new5')
ggsave(paste0(name, '.pdf'), plot = p, width = 35, height = 25, units = "cm", dpi = 300)
ggsave(paste0(name, '.png'), plot = p, width = 35, height = 25, units = "cm", dpi = 300)



# Map of the bias that works
# Creating the dataframe
# For pr
bias_tmp <- GC_hellinger_projections$pr - reference_list$future$pr[[1]]

# Ensure you are removing NA values in the computation
rm(avg_bias)
avg_bias <- round(mean(abs(bias_tmp), na.rm = TRUE), 2)

# Settings the limits
limit <- c(-10, 10)

# Settings the limits

v_limits <- as.numeric(format(seq(min(limit), max(limit), len=5), digits = 3))
bias_tmp[bias_tmp < min(limit)] <- min(limit)
bias_tmp[bias_tmp > max(limit)] <- max(limit)

test_df <- melt(bias_tmp, c("lon", "lat"), value.name = "Bias")

p <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias),)+
  ggtitle(paste0('GraphCut Hellinger'))+
  labs(subtitle = paste0('Average Bias : ', avg_bias))+
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),
                       breaks = v_limits,
                       limits = limit)+
  # guide = guide_colourbar(
  #   barwidth = 0.4,
  #   ticks.colour = 'black',
  #   ticks.linewidth = 1,
  #   frame.colour = 'black',
  #   draw.ulim = TRUE,
  #   draw.llim = TRUE),)
  borders("world2", colour = 'black', lwd = 0.12) +
  scale_x_continuous(, expand = c(0, 0)) +
  scale_y_continuous(, expand = c(0,0))+
  theme(legend.position = 'bottom')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())+
  xlab('Longitude')+
  ylab('Latitude') +
  labs(fill='[mm/day]')+
  theme_bw()+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=12))+ #change legend text font size
  theme(plot.title = element_text(size=24),
        plot.subtitle = element_text(size = 20,hjust=0.5),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),)+
  easy_center_title()
p
mean(abs(bias_tmp))

name <- paste0('figure/GC_Hellinger_new_bias_pr_26models_new5')
ggsave(paste0(name, '.pdf'), plot = p, width = 35, height = 25, units = "cm", dpi = 300)
ggsave(paste0(name, '.png'), plot = p, width = 35, height = 25, units = "cm", dpi = 300)





h_dist_map <- array(NA, dim = dim(GC_result_hellinger_test_new_era5_5$label_attribution))

for(j in seq_along(variables)){
  for(l in 0:(length(model_names))){
    islabel <- which(GC_result_hellinger_test_new_era5_5$label_attribution == l)
    h_dist_map[islabel] <- h_dist[,,(l)][islabel]
  }
}

test_df <- melt(h_dist_map, c("lon", "lat"), value.name = "Bias")

p5 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias))+
  labs(subtitle = 'Projection period : 1977 - 1999')+
  ggtitle(paste0('GraphCut Hellinger', ': Mean Hellinger distance = ', round(mean(h_dist_map), 2)))+
  scale_fill_gradient(low = "white", high = "#015a8c", limits = c(0.1, 0.70), oob = scales::squish)+
  borders("world2", colour = 'black', lwd = 0.12) +
  scale_x_continuous(, expand = c(0, 0)) +
  scale_y_continuous(, expand = c(0,0))+
  theme(legend.position = 'bottom')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())+
  xlab('Longitude')+
  ylab('Latitude') +
  labs(fill='Hellinger \nDistance')+
  theme_bw()+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=16), #change legend title font sizen
        legend.text = element_text(size=12))+ #change legend text font size
  theme(plot.title = element_text(size=24),
        plot.subtitle = element_text(size = 20,hjust=0.5),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),)+
  easy_center_title()
p5

name <- paste0('figure/H_dist_future_GC_hellinger_new_26models_present5')
ggsave(paste0(name, '.pdf'), plot = p5, width = 35, height = 25, units = "cm", dpi = 300)
ggsave(paste0(name, '.png'), plot = p5, width = 35, height = 25, units = "cm", dpi = 300)




