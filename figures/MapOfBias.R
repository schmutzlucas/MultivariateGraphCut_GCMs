gc_bias_var <- bias_var(
  var_future = list(models_matrix$future[,,, 2], models_matrix$future[,,, 1]),
  ref_future_nrm = list(reference_matrix_nrm$future[,,1,2], reference_matrix_nrm$future[,,1,1]),
  var_future_nrm = list(models_matrix_nrm$future[,,,2], models_matrix_nrm$future[,,,1]),
  labeling = GC_result_hellinger$label_attribution
)

# Map of the bias that works
# Creating the dataframe
test_df <- melt(gc_bias_var$bias[[1]], c("lon", "lat"), value.name = "Bias")
test_df$Bias <- test_df$Bias * meta_normalization_ToP$tas_future$sd

# Settings the limits
limit <- ceiling(max(abs(c(min(test_df$Bias), max(test_df$Bias)))))
limit <- round_any(limit, 10, floor)
limit <- 8
limits <- c(-8, 8)
v_limits <- as.numeric(format(seq(-limit, limit, len=6), digits = 3))
test_df$Bias[test_df$Bias < -limit] <- -limit
test_df$Bias[test_df$Bias > limit] <- limit

p <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat, fill=Bias),)+
  ggtitle('Biases in GC hybrid for temperatures')+
  scale_fill_gradientn(colours = rev((brewer.rdbu(100))),
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
  labs(fill='[K]')+
  theme_bw()+
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.7, 'cm'), #change legend key height
        legend.key.width = unit(0.2, 'cm'), #change legend key width
        legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=7))+ #change legend text font size
  theme(plot.title = element_text(size=16),
        axis.text=element_text(size=7),
        axis.title=element_text(size=9),)+
  easy_center_title()
p