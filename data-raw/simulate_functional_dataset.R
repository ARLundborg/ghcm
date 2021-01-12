library(pracma)
library(refund)

evaluate_on_grid <- function(list_of_functions, grid) {
  #' Evaluates the list of functions on the given grid and return the output in a matrix
  #' with each row representing a function and each column representing a grid point

  list_of_evaluated_functions <- lapply(list_of_functions, function(f) {f(grid)})
  k <- length(list_of_evaluated_functions)
  matrix(unlist(list_of_evaluated_functions, use.names = F), ncol = k)
}

brownian_noise <- function(N, grid, sigma=1) {
  #' Returns N Brownian motions with mean zero and variance given by sigma^2.
  #' Each row of the outputted matrix corresponds to a single Brownian motion
  #'  and each column corresponds to a grid point.
  diff_grid <- diff(grid)
  t(replicate(N, cumsum(c(0, rnorm(length(diff_grid), 0, sd = sigma*sqrt(diff_grid))))))
}

linear_operator <- function(grid, X, beta) {
  #' Evaluate a linear operator of the form \int \beta(s, t) X(t) \, \mathrm{d}t given X on a grid.
  #' Each row of the outputted matrix corresponds to the linear operator applied to each row of X and
  #' each column corresponds to a grid point. It is assumed that X is observed on grid.

  N <- dim(X)[1]
  beta_vals <- outer(grid, grid, FUN = beta)
  t(sapply(1:N, function(i){
    sapply(1:length(grid), function(j) {
      trapz(grid, X[i, ]*beta_vals[j, ])
    })
  }))
}

linear_operator_1d <- function(grid, X, alpha) {
  #' Evaluate a linear operator of the form \int \alpha(t) X(t) \, \mathrm{d}t given X on a grid.
  #' The output is a vector where the i'th entry corresponds to the linear operator applied to the i'th row of X.
  #' It is assumed that X is observed on grid.

  N <- dim(X)[1]
  alpha_vals <- alpha(grid)
  sapply(1:N, function(i) trapz(grid, X[i, ] * alpha_vals))
}


set.seed(011221)

N <- 500
grid <- seq(0, 1, length.out=101)
beta_Z <- function(s, t) sin(2*pi*s*t)+cos(2*pi*(1-s*t)) + s^2+t^2
sigma_Z <- 0.5
beta_W <- function(s, t) sin(6*pi*s*t)*sqrt(s*t)
sigma_W <- 0.15
alpha <- function(t) t^2
sigma_1 <- 0.20
sigma_2 <- 0.20

X <- brownian_noise(N, grid)
Z <- linear_operator(grid, X, beta_Z) + brownian_noise(N, grid, sigma_Z)
W <- linear_operator(grid, Z, beta_W) + brownian_noise(N, grid, sigma_W)
Y_1 <- linear_operator_1d(grid, Z, alpha) + rnorm(N, sd = sigma_1)
Y_2 <- linear_operator_1d(grid, Z, alpha) + rnorm(N, sd = sigma_2)


ghcm_sim_data <- data.frame(Y_1=Y_1, Y_2=Y_2)
ghcm_sim_data$X <- X
ghcm_sim_data$Z <- Z
ghcm_sim_data$W <- W

usethis::use_data(ghcm_sim_data, overwrite = TRUE)
