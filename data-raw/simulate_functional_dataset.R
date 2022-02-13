library(pracma)
library(refund)

evaluate_on_grid <- function(list_of_functions, grid) {
  #' Evaluates the list of functions on the given grid and return the output in a matrix
  #' with each row representing a function and each column representing a grid point

  list_of_evaluated_functions <- lapply(list_of_functions, function(f) {f(grid)})
  k <- length(list_of_evaluated_functions)
  matrix(unlist(list_of_evaluated_functions, use.names = F), ncol = k)
}

brownian_noise <- function(n, grid, sigma=1) {
  #' Returns n Brownian motions with mean zero and variance given by sigma^2.
  #' Each row of the outputted matrix corresponds to a single Brownian motion
  #'  and each column corresponds to a grid point.
  diff_grid <- diff(grid)
  t(replicate(n, cumsum(c(0, rnorm(length(diff_grid), 0, sd = sigma*sqrt(diff_grid))))))
}

linear_operator <- function(grid, X, beta) {
  #' Evaluate a linear operator of the form \int \beta(s, t) X(t) \, \mathrm{d}t given X on a grid.
  #' Each row of the outputted matrix corresponds to the linear operator applied to each row of X and
  #' each column corresponds to a grid point. It is assumed that X is observed on grid.

  n <- dim(X)[1]
  beta_vals <- outer(grid, grid, FUN = beta)
  t(sapply(1:n, function(i){
    sapply(1:length(grid), function(j) {
      trapz(grid, X[i, ]*beta_vals[j, ])
    })
  }))
}

linear_operator_1d <- function(grid, X, beta) {
  #' Evaluate a linear operator of the form \int \beta(t) X(t) \, \mathrm{d}t given X on a grid.
  #' The output is a vector where the i'th entry corresponds to the linear operator applied to the i'th row of X.
  #' It is assumed that X is observed on grid.

  n <- dim(X)[1]
  beta_vals <- beta(grid)
  sapply(1:n, function(i) trapz(grid, X[i, ] * beta_vals))
}

sub_sample_regular_functional_data <- function(X, grid, lambda) {
  #' Given a \code{grid} of size k and a n x k matrix, \code{X}, of functional
  #' observations, this function subsamples each row of X. A Poisson random
  #' variable with mean \code{lambda} is generated for each row determining the
  #' number of grid points to be subsampled (although at least 4 are always
  #' sampled).
  #'
  #' The data is returned as a "melted" data.frame with 3 columns; \code{.obs}
  #' indicating which curve the row represents, \code{.index} giving the
  #' function argument and \code{.value} containing the function value.
  #'
  #'
  #'

  n <- dim(X)[1]
  number_of_grid_points <- length(grid)
  df_list <- lapply(1:n, function(i) {
    idx <- sample(1:number_of_grid_points, max(rpois(1, lambda), 4))
    grid_points <- grid[idx]
    values <- X[i, idx]
    data.frame(.obs=i, .index=grid_points, .value=values)
  })
  return(Reduce(rbind.data.frame, df_list))
}


set.seed(011221)

n <- 500
grid <- seq(0, 1, length.out=101)
beta_Z <- function(s, t) sin(2*pi*s*t)+cos(2*pi*(1-s*t)) + s^2+t^2
sigma_Z <- 0.5
beta_W <- function(s, t) sin(6*pi*s*t)*sqrt(s*t)
sigma_W <- 0.15
alpha <- function(t) t^2
sigma_1 <- 0.20
sigma_2 <- 0.20
lambda_X <- 10
lambda_W <- 15

X <- brownian_noise(n, grid)
Z <- linear_operator(grid, X, beta_Z) + brownian_noise(n, grid, sigma_Z)
W <- linear_operator(grid, Z, beta_W) + brownian_noise(n, grid, sigma_W)
Y_1 <- linear_operator_1d(grid, Z, alpha) + rnorm(n, sd = sigma_1)
Y_2 <- linear_operator_1d(grid, Z, alpha) + rnorm(n, sd = sigma_2)


ghcm_sim_data <- data.frame(Y_1=Y_1, Y_2=Y_2)
ghcm_sim_data$X <- X
ghcm_sim_data$Z <- Z
ghcm_sim_data$W <- W

usethis::use_data(ghcm_sim_data, overwrite = TRUE)

X_df <- sub_sample_regular_functional_data(X, grid, lambda_X)
W_df <- sub_sample_regular_functional_data(W, grid, lambda_W)

ghcm_sim_data_irregular <- list(Y_1=Y_1, Y_2=Y_2, Z=Z, X=X_df, W=W_df)

usethis::use_data(ghcm_sim_data_irregular, overwrite = TRUE)


