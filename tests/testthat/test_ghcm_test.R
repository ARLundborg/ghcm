context("ghcm_test")

library(refund)
set.seed(13720)

N <- 200
grid_X <- seq(0, 1, length.out = 50)
grid_Y <- c(seq(0, 0.5, length.out = 45), 0.6, 0.7, 0.8, 0.9, 1)
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
  tst <- ghcm_test(resid_X_on_Z_multivariate, resid_Y_on_Z_scalar, X_grid = NA)
  expect_is(tst, "ghcm")
})

test_that("ghcm_runs_with_functional_and_scalar_inputs", {
  tst <- ghcm_test(resid_X_on_Z_functional, resid_Y_on_Z_scalar,
                   X_grid = grid_X)
  expect_is(tst, "ghcm")
})

test_that("ghcm_runs_with_two_functional_inputs", {
  tst <- ghcm_test(resid_X_on_Z_functional, resid_Y_on_Z_functional,
                   X_grid = grid_X, Y_grid = grid_Y)
  expect_is(tst, "ghcm")
})

test_that("ghcm_fails_when_sample_sizes_differ", {
  resid_X_on_Z_functional_half <- resid_X_on_Z_functional[1:(N/2), ]

  expect_error({
    ghcm_test(resid_X_on_Z_functional_half, resid_Y_on_Z_functional,
              X_grid = grid_X, Y_grid = grid_Y)
  })
})

test_that("ghcm_fails_when_giving_bad_fpca_string", {
  expect_error({
    ghcm_test(resid_X_on_Z_functional, resid_Y_on_Z_functional,
              X_grid = grid_X, Y_grid = grid_Y,
              fpca_method = "not_an_fpca_method")
  })
})

test_that("ghcm_warns_about_X_non_equidistant_grid", {
  expect_warning({
    ghcm_test(resid_X_on_Z_functional, resid_Y_on_Z_scalar)
  }, "equidistant")
  tst <- ghcm_test(resid_X_on_Z_functional, resid_Y_on_Z_scalar,
                   X_grid = grid_X)
  expect_is(tst, "ghcm")
})

test_that("ghcm_warns_about_Y_non_equidistant_grid", {
  expect_warning({
    ghcm_test(resid_X_on_Z_scalar, resid_Y_on_Z_functional)
  }, "equidistant")
  tst <- ghcm_test(resid_X_on_Z_scalar, resid_Y_on_Z_functional,
                   Y_grid = grid_Y)
  expect_is(tst, "ghcm")
})
