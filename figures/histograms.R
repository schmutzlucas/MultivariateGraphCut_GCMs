
# Calculate the maximum value for the x-axis limit
max_value <- max(GC_result_hellinger$label_attribution)

# Create a sequence of breaks at every half unit starting from -0.5
breaks <- seq(-0.5, max_value + 0.5, by = 1)

# Generate the histogram with specified breaks
hist(GC_result_hellinger$label_attribution, breaks = breaks, xlim = c(0, max_value), xaxt = 'n', freq = FALSE)

# Add custom x-axis ticks centered on the bars
axis(1, at = seq(0, max_value, by = 1), labels = seq(0, max_value, by = 1))


# Calculate the maximum value for the x-axis limit
max_value <- max(GC_result_hellinger$label_attribution)

# Create a sequence of breaks at every half unit starting from -0.5
breaks <- seq(-0.5, max_value + 0.5, by = 1)

# Generate the histogram with specified breaks
hist(GC_result_hellinger$label_attribution, breaks = breaks, xlim = c(0, max_value), xaxt = 'n', freq = FALSE, main = "Frequency of label attribution", xlab = "Model index", ylab = "")

# Add custom x-axis ticks centered on the bars
axis(1, at = seq(0, max_value, by = 1), labels = seq(0, max_value, by = 1))



# Calculate the maximum value for the x-axis limit
max_value <- max(GC_result_hellinger$label_attribution)

# Create a sequence of breaks at every half unit starting from -0.5
breaks <- seq(-0.5, max_value + 0.5, by = 1)

# Generate the histogram with specified breaks
hist(GC_result_hellinger$label_attribution, breaks = breaks, xlim = c(0, max_value), xaxt = 'n', freq = FALSE, main = "Density of label attribution", xlab = "Model index", ylab = "")

# Add custom x-axis ticks centered on the bars
axis(1, at = seq(0, max_value, by = 1), labels = seq(0, max_value, by = 1))





