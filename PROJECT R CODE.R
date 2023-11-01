### Start measuring time
start_time <- Sys.time()

####### Parameters we are going to use.

### Initial stock price
S0 <- 100
### Strike price
K <- 100
### Time to maturity
T <- 1 
### Risk-free rate
r <- 0.05 
### Volatility
sigma <- 0.2
### Barrier level
B <- 105 
### Number of simulations
N <- 50000 
### Number of time steps
M <- 100
### The time increment
dt <- T / M
  
# Function to simulate stock price paths
set.seed(42)
simulate_paths <- function(N, M, S0, r, sigma, dt) {
  paths <- matrix(0, nrow = N, ncol = M + 1)
  paths[, 1] <- S0
  
  for (i in 1:N) {
    for (j in 2:(M + 1)) {
      paths[i, j] <- paths[i, j - 1] * exp((r - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * rnorm(1))
    }
  }
  
  return(paths)
}

# Function to simulate barrier option paths
simulate_barrier_paths <- function(N, M, S0, r, sigma, dt, B, type) {
  paths <- simulate_paths(N, M, S0, r, sigma, dt)
  crossed_barrier <- apply(paths > B, 2, any)
  
  if (type == "in") {
    paths[1, !crossed_barrier] <- S0
  } else if (type == "out") {
    paths[1, crossed_barrier] <- S0
  }
  
  return(paths)
}

# Simulate stock price paths for each scenario
down_and_in_paths <- simulate_barrier_paths(1, M, S0, r, sigma, dt, B, "in")
down_and_out_paths <- simulate_barrier_paths(1, M, S0, r, sigma, dt, B, "out")
up_and_in_paths <- simulate_barrier_paths(1, M, S0, r, sigma, dt, B, "in")
up_and_out_paths <- simulate_barrier_paths(1, M, S0, r, sigma, dt, B, "out")

# Visualization with different colors
#par(mfrow = c(2, 2))
plot(down_and_in_paths[1, ], type = "l", col = "blue", xlab = "Time Steps", ylab = "Stock Price", main = "Down-and-In Option")
abline(h = B, col = "red", lty = 2)
plot(down_and_out_paths[1, ], type = "l", col = "green", xlab = "Time Steps", ylab = "Stock Price", main = "Down-and-Out Option")
abline(h = B, col = "red", lty = 2)
plot(up_and_in_paths[1, ], type = "l", col = "red", xlab = "Time Steps", ylab = "Stock Price", main = "Up-and-In Option")
abline(h = B, col = "blue", lty = 2)
plot(up_and_out_paths[1, ], type = "l", col = "purple", xlab = "Time Steps", ylab = "Stock Price", main = "Up-and-Out Option")
abline(h = B, col = "blue", lty = 2)

#####################################################################################

###### Simulation using the manipulating log-normal model
#set.seed(42)
paths <- matrix(0, nrow=N, ncol=M+1)
paths[, 1] <- S0
for (i in 1:N) {
  for (j in 2:(M+1)) {
    paths[i, j] <- paths[i, j-1] * exp((r - 0.5 * sigma^2) * dt +
                                         sigma * sqrt(dt) * rnorm(1))
  }
}

### Visualization of Simulated Stock Price Paths
### Number of paths to plot for visualization
num_paths_to_plot <- 100  

### Colors for paths
colors <- ifelse(apply(paths[1:num_paths_to_plot, 2:(M+1)], 1, function(x) any(x > B)), "red", "blue")

matplot(t(paths[1:num_paths_to_plot, ]), type="l", lty=1, xlab="Time Steps", ylab="Stock Price", main="Simulated Stock Price Paths", col=colors)
abline(h=B, col="black", lty=2) 
legend("topright", legend=c("Crossed Barrier", "Did Not Cross"), fill=c("red", "blue"))


### Calculate payoffs
payoffs <- pmax(paths[, M+1] - K, 0)

### WE identify paths that crossed the barrier level
crossed_barrier <- apply(paths[, 2:(M+1)] > B, 1, any)
num_crossed <- sum(crossed_barrier)
num_not_crossed <- N - num_crossed

cat("Number of simulations that crossed the barrier:", num_crossed, "\n")
cat("Number of simulations that did not cross the barrier:", num_not_crossed, "\n")

########################################################################

## WE ARE GOING TO CONSIDER LOOK AT FIVE DIFFERENT TYPES FINITE BASIS FUNCTIONS 
## THEN WE CHOOSE WHICH ONE PERFORM BETTER FOR OUR PROBLEM.

################# THIN-PLATE SPLINES FINITE BASIS FUNCTION #####################
### Load necessary libraries
library(fields)
library(mgcv)

### Calculate payoffs for European up-and-in call option
### Indicator function for paths that crossed the barrier
I_B <- ifelse(crossed_barrier, 1, 0)

### Payoff for each path at maturity
payoff_at_maturity <- pmax(paths[, M+1] - K, 0)

### Compute the option price using the formula
option_price_direct <- exp(-r * T) * mean(I_B * payoff_at_maturity)

### Select paths that crossed the barrier
X_crossed <- paths[crossed_barrier, M]
payoffs_crossed <- payoffs[crossed_barrier]

### Least Squares Regression using Thin-Plate Splines
model <- gam(payoffs_crossed ~ s(X_crossed, bs="tp"))

### Option Price Estimation using the regression model
estimated_payoffs <- predict(model, newdata=data.frame(X_crossed=X_crossed))

### Multiply the estimated payoffs by the indicator function I_B
final_payoffs <- I_B[crossed_barrier] * estimated_payoffs

### Compute the option price using the direct formula
option_price <- exp(-r * T) * mean(final_payoffs)

### Residual Analysis
residuals <- residuals(model)
hist(residuals, breaks=50, main="Histogram of Residuals", xlab="Residuals")
plot(residuals, type="p", main="Residual Plot", xlab="Index", ylab="Residuals")
abline(h=0, col="red")

### Model Performance Metrics
MAE <- mean(abs(estimated_payoffs - payoffs_crossed))
RMSE <- sqrt(mean((estimated_payoffs - payoffs_crossed)^2))
### Compute R-squared manually
observed <- payoffs_crossed
predicted <- estimated_payoffs
RSS <- sum((observed - predicted)^2)
TSS <- sum((observed - mean(observed))^2)
R_squared <- 1 - (RSS / TSS)

### End measuring time
end_time <- Sys.time()

### Calculate the elapsed time
elapsed_time <- end_time - start_time

### Print the results
cat("Elapsed time:", elapsed_time, "\n")
print(paste("Directly Computed Option Price:", option_price_direct))
print(paste("Estimated Option Price:", option_price))
print(paste("Mean Absolute Error (MAE):", MAE))
print(paste("Root Mean Square Error (RMSE):", RMSE))
print(paste("R-squared (R_squared):", R_squared))


######################################################################
################# B-SPLINES FINITE BASIS FUNCTION #####################

### Load necessary libraries
library(splines)

### Basis functions using B-splines
basis_functions <- function(x) {
  bs(x, degree=3, knots=c(75, 85, 95, 105, 115))  
}
X <- basis_functions(paths[, M])

### Calculate payoffs for European up-and-in call option
### Indicator function for paths that crossed the barrier
I_B <- ifelse(crossed_barrier, 1, 0)

### Payoff for each path at maturity
payoff_at_maturity <- pmax(paths[, M+1] - K, 0)

### Compute the option price using the formula
option_price_direct <- exp(-r * T) * mean(I_B * payoff_at_maturity)

### Least Squares Regression using built-in lm function
model <- lm(payoff_at_maturity ~ X - 1)

### Option Price Estimation using the regression model
X_crossed <- X[crossed_barrier, ]
estimated_payoffs <- predict(model, newdata=data.frame(X_crossed))

### Multiply the estimated payoffs by the indicator function I_B
final_payoffs <- I_B[crossed_barrier] * estimated_payoffs

### Compute the option price using the direct formula
option_price <- exp(-r * T) * mean(final_payoffs)

### Residual Analysis
### Plot residuals vs. fitted values
plot(model, which = 1)
residuals <- residuals(model)
hist(residuals, breaks=50, main="Histogram of Residuals", xlab="Residuals")
plot(residuals, type="p", main="Residual Plot", xlab="Index", ylab="Residuals")
abline(h=0, col="red")

### Model Performance Metrics
MAE <- mean(abs(estimated_payoffs - payoff_at_maturity))
RMSE <- sqrt(mean((estimated_payoffs - payoff_at_maturity)^2))
R_squared <- summary(model)$r.squared

### End measuring time
end_time <- Sys.time()

### Calculate the elapsed time
elapsed_time <- end_time - start_time

### Print the results
cat("Elapsed time:", elapsed_time, "\n")
print(paste("Directly Computed Option Price:", option_price_direct))
print(paste("Estimated Option Price:", option_price))
print(paste("Mean Absolute Error (MAE):", MAE))
print(paste("Root Mean Square Error (RMSE):", RMSE))
print(paste("R-squared (R_squared):", R_squared))

#############################################################################

########### SPARSE PIECEWISE LINEAR Functions FINITE BASIS FUNCTION #############

### Sparse Piecewise Linear Functions as the basis function
basis_functions <- function(x) {
  breakpoints <- c(80, 90, 100, 110)  
  as.numeric(cut(x, breaks=breakpoints, labels=FALSE, include.lowest=TRUE))
}
X <- matrix(basis_functions(paths[, M]), ncol=1)

### Calculate payoffs for European up-and-in call option
### Indicator function for paths that crossed the barrier
I_B <- ifelse(crossed_barrier, 1, 0)

### Payoff for each path at maturity
payoff_at_maturity <- pmax(paths[, M+1] - K, 0)

### Compute the option price using the formula
option_price_direct <- exp(-r * T) * mean(I_B * payoff_at_maturity)

### Least Squares Regression
### Using built-in lm function for regression
model <- lm(payoff_at_maturity ~ X - 1)

### Option Price Estimation using the regression model
X_crossed <- X[crossed_barrier, ]
estimated_payoffs <- predict(model, newdata=data.frame(X_crossed=X_crossed))

### Multiply the estimated payoffs by the indicator function I_B
final_payoffs <- I_B[crossed_barrier] * estimated_payoffs

### Compute the option price using the direct formula
option_price <- exp(-r * T) * mean(final_payoffs)

### Residual Analysis
residuals <- residuals(model)
hist(residuals, breaks=50, main="Histogram of Residuals", xlab="Residuals")
plot(residuals, type="p", main="Residual Plot", xlab="Index", ylab="Residuals")
abline(h=0, col="red")

### Model Performance Metrics
MAE <- mean(abs(estimated_payoffs - payoff_at_maturity))
RMSE <- sqrt(mean((estimated_payoffs - payoff_at_maturity)^2))
R_squared <- summary(model)$r.squared

### End measuring time
end_time <- Sys.time()

### Calculate the elapsed time
elapsed_time <- end_time - start_time

### Print the results
cat("Elapsed time:", elapsed_time, "\n")
print(paste("Directly Computed Option Price:", option_price_direct))
print(paste("Estimated Option Price:", option_price))
print(paste("Mean Absolute Error (MAE):", MAE))
print(paste("Root Mean Square Error (RMSE):", RMSE))
print(paste("R-squared (R_squared):", R_squared))

###################################################################

############## SPARSE POLYNOMIAL FINITE BASIS FUNCTION #######################

### Basis functions using Sparse Polynomial for Only constant, linear, and cubic terms
basis_functions <- function(x) {
  cbind(1, x, x^3)  
}
X <- basis_functions(paths[, M])

### Calculate payoffs for European up-and-in call option
### Indicator function for paths that crossed the barrier
I_B <- ifelse(crossed_barrier, 1, 0)

### Payoff for each path at maturity
payoff_at_maturity <- pmax(paths[, M+1] - K, 0)

### Compute the option price using the formula
option_price_direct <- exp(-r * T) * mean(I_B * payoff_at_maturity)

### Least Squares Regression
### Using built-in lm function for regression
model <- lm(payoff_at_maturity ~ X - 1)

### Option Price Estimation using the regression model
X_crossed <- X[crossed_barrier, ]
estimated_payoffs <- predict(model, newdata=data.frame(X_crossed=X_crossed))

### Multiply the estimated payoffs by the indicator function I_B
final_payoffs <- I_B[crossed_barrier] * estimated_payoffs

### Compute the option price using the direct formula
option_price <- exp(-r * T) * mean(final_payoffs)

### Residual Analysis
residuals <- residuals(model)
hist(residuals, breaks=50, main="Histogram of Residuals", xlab="Residuals")
plot(residuals, type="p", main="Residual Plot", xlab="Index", ylab="Residuals")
abline(h=0, col="red")

### Model Performance Metrics
MAE <- mean(abs(estimated_payoffs - payoff_at_maturity))
RMSE <- sqrt(mean((estimated_payoffs - payoff_at_maturity)^2))
R_squared <- summary(model)$r.squared

### End measuring time
end_time <- Sys.time()

### Calculate the elapsed time
elapsed_time <- end_time - start_time

### Print the results
cat("Elapsed time:", elapsed_time, "\n")
print(paste("Directly Computed Option Price:", option_price_direct))
print(paste("Estimated Option Price:", option_price))
print(paste("Mean Absolute Error (MAE):", MAE))
print(paste("Root Mean Square Error (RMSE):", RMSE))
print(paste("R-squared (R^2):", R_squared))

########################################################################

############## POLYNOMIAL FINITE BASIS FUNCTION #######################

### Basis functions
### Using a simpler set of basis functions
basis_functions <- function(x) {
  cbind(1, x, x^2)  
}
X <- basis_functions(paths[, M])

### Least Squares Regression
### Using built-in lm function for regression
model <- lm(payoff_at_maturity ~ X - 1)

### Option Price Estimation using the regression model
X_crossed <- X[crossed_barrier, ]
estimated_payoffs <- predict(model, newdata=data.frame(X_crossed=X_crossed))

### Multiply the estimated payoffs by the indicator function I_B
final_payoffs <- I_B[crossed_barrier] * estimated_payoffs

### Compute the option price using the direct formula
option_price <- exp(-r * T) * mean(final_payoffs)

### Residual Analysis
residuals <- residuals(model)
hist(residuals, breaks=50, main="Histogram of Residuals", xlab="Residuals")
plot(residuals, type="p", main="Residual Plot", xlab="Index", ylab="Residuals")
abline(h=0, col="red")

### Model Performance Metrics
MAE <- mean(abs(estimated_payoffs - payoff_at_maturity))
RMSE <- sqrt(mean((estimated_payoffs - payoff_at_maturity)^2))
R_squared <- summary(model)$r.squared

### End measuring time
end_time <- Sys.time()

### Calculate the elapsed time
elapsed_time <- end_time - start_time

### Print the results
cat("Elapsed time:", elapsed_time, "\n")
print(paste("Directly Computed Option Price:", option_price_direct))
print(paste("Estimated Option Price:", option_price))
print(paste("Mean Absolute Error (MAE):", MAE))
print(paste("Root Mean Square Error (RMSE):", RMSE))
print(paste("R-squared (R^2):", R_squared))



#################################################################
### BAR GRAPH OF THE RESULTS OBTAIN FROM THE ABOVE MODELS USING R^2 VAULES

### Load necessary libraries
library(ggplot2)
library(dplyr)

### Create a data frame from the table
data <- data.frame(
  Finite_Basis_Function = c("Polynomial", "B-splines", "Sparse polynomial", "Sparse piecewise linear functions", "Thin-plate splines"),
  R_squared = c(0.9709, 0.9909, 0.9635, 0.4614, 0.9846)
)

### Create the plot
ggplot(data, aes(x = Finite_Basis_Function, y = R_squared, fill = Finite_Basis_Function)) +
  geom_bar(stat = "identity") +
  labs(title = "Models Performance based on R-squared", y = "R-squared (R^2)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Sparse polynomial" = "red", "Polynomial" = "blue", "B-splines" = "green", "Thin-plate splines" = "purple", "Sparse piecewise linear functions" = "orange"))
