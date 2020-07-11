sim_brownian_motion <- function(n, grid, sigma=1) {
  #' Simulate a Brownian Motion
  #'
  #' \code{sim_brownian_motion} returns \eqn{n} simulated Brownian motion sample
  #' paths evaluated on a grid.
  #'
  #' @param n Positive integer, the number of sample paths to simulate.
  #' @param grid Numeric vector with consecutive values in the unit interval,
  #'   the grid points to evalute the sample paths.
  #' @param sigma Optionally, a positive number specifying the square root of
  #'   the variance of the Brownian motions to simulate. Defaults to 1.
  #'
  #' @return Returns an \eqn{n} x \code{len(grid)} matrix where each row
  #' is a sample path evaluated at the grid points.
  #'
  #' @examples
  #' grid <- seq(0, 1, length.out = 100)
  #' n <- 1000
  #'
  #' X <- ghcm:::sim_brownian_motion(n, grid)
  #'
  #' sigma <- 2
  #'
  #' Y <- ghcm:::sim_brownian_motion(n, grid, sigma)
  #'

  diff_grid <- diff(grid)
  t(replicate(n, cumsum(c(0, stats::rnorm(length(diff_grid), 0,
                                   sd = sigma * sqrt(diff_grid))))))
}
