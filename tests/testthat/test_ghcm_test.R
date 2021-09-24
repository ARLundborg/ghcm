context("ghcm_test")

library(refund)
set.seed(13720)

N <- 200
grid_X <- seq(0, 1, length.out = 50)
grid_Y <- c(seq(0, 2, length.out = 25))
resid_X_on_Z_scalar <- rnorm(N)
resid_X_on_Z_multivariate <- matrix(rnorm(2 * N), ncol = 2)
resid_X_on_Z_functional <- matrix(rnorm(length(grid_X) * N),
                                  ncol = length(grid_X))
resid_Y_on_Z_scalar <- rnorm(N)
resid_Y_on_Z_functional <- matrix(rnorm(length(grid_Y) * N),
                                  ncol = length(grid_Y))

test_that("ghcm_runs_with_two_scalar_inputs", {
  tst <- ghcm_test(resid_X_on_Z_scalar, resid_Y_on_Z_scalar)
  expect_is(tst, "ghcm")
})

test_that("ghcm_runs_with_multivariate_and_scalar_inputs", {
  tst <- ghcm_test(resid_X_on_Z_multivariate, resid_Y_on_Z_scalar)
  expect_is(tst, "ghcm")
})

test_that("ghcm_runs_with_functional_and_scalar_inputs", {
  tst <- ghcm_test(resid_X_on_Z_functional, resid_Y_on_Z_scalar)
  expect_is(tst, "ghcm")
})

test_that("ghcm_runs_with_two_functional_inputs", {
  tst <- ghcm_test(resid_X_on_Z_functional, resid_Y_on_Z_functional)
  expect_is(tst, "ghcm")
})

test_that("ghcm_fails_when_sample_sizes_differ", {
  resid_X_on_Z_functional_half <- resid_X_on_Z_functional[1:(N/2), ]

  expect_error({
    ghcm_test(resid_X_on_Z_functional_half, resid_Y_on_Z_functional)
  })
})
