kappa <- 1.5
# Define a function to calculate the average cancer rate based on a parameter `lambda`
avg_cancer_rate <- function(lambda) {
  # Compute the difference between two exponential survival functions,
  # representing the probability of surviving without cancer up to two specific ages (40.9 and 51.4 years, converted to days).
  exp(-((40.9225 * 365/lambda)^kappa)) - exp(-((51.4225 * 365/lambda)^kappa))
}

# Define an error function that measures the squared difference between the average cancer rate and a target value (0.1).
err_func <- function(lambda) {
  (avg_cancer_rate(lambda) - 0.1)^2
}

# Use the `optimise` function to find the value of `lambda` that minimizes the error function within a range (0 to 1,000,000).
# The result will provide the optimal lambda where the average cancer rate is closest to 0.1.
optimise(err_func, c(0, 1000000), tol = 1e-10, maximum = FALSE)

# Create a sequence of lambda values to evaluate the functions across a range.
lambda <- seq(0, 100000, length.out = 10000)

# Plot the average cancer rate as a function of lambda, showing how it changes across the specified range.
plot(lambda, avg_cancer_rate(lambda), type = "l", xlab = "lambda", ylab = "Average cancer rate")

# Plot the error function as a function of lambda to visualize the error across the range and identify the optimal value.
plot(lambda, err_func(lambda), type = "l", xlab = "lambda", ylab = "Error")
