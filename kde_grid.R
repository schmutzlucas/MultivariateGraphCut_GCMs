library(ks)


# Number of points at which to evaluate the pdf
num_pdf_points <- 100

# Initialize the output array with dimensions 360 x 181 x num_pdf_points
pdf_output <- array(0, dim = c(360, 181, num_pdf_points))

# Iterate through the grid points
for (lon in 1:360) {
  for (lat in 1:181) {
    # Extract the time series data for the current grid point
    time_series <- test[lon, lat, ]

    # Compute the kernel density estimate for the time series
    kde <- kde(time_series)

    # Define the points at which to evaluate the pdf
    x_points <- seq(min(time_series), max(time_series), length.out = num_pdf_points)

    # Compute the pdf values at the defined points
    pdf_values <- predict(kde, x = x_points)

    # Store the pdf values in the output array
    pdf_output[lon, lat, ] <- pdf_values
  }
}

# The 'pdf_output' array now contains the desired output
