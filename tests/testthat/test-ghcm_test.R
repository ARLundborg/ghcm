context("ghcm_test")

library(refund)
set.seed(13720)

n <- 200
grid_X <- seq(0, 1, length.out = 50)
grid_Y <- c(seq(0, 2, length.out = 25))

resid_X_on_Z_scalar <- rnorm(n)
resid_X_on_Z_multivariate <- matrix(rnorm(2 * n), ncol = 2)
resid_X_on_Z_irregular <- Reduce(rbind.data.frame, lapply(1:n, function(i) {
  k <- max(rpois(1, 10), 4)
  data.frame(.obs=i, .index=runif(k), .value=rnorm(k))
}))

resid_Y_on_Z_scalar <- rnorm(n)
resid_Y_on_Z_functional <- matrix(rnorm(length(grid_Y) * n),
                                  ncol = length(grid_Y))
resid_Y_on_Z_irregular <- Reduce(rbind.data.frame, lapply(1:n, function(i) {
  k <- max(rpois(1, 6), 4)
  data.frame(.obs=i, .index=runif(k, 0, 2), .value=rnorm(k))
}))

test_that("ghcm_runs_with_two_irregular_inputs", {
  tst <- ghcm_test(resid_X_on_Z_irregular, resid_Y_on_Z_irregular, c(0,1), c(0,2))
  expect_is(tst, "ghcm")
})

test_that("ghcm_runs_with_one_regular_and_one_irregular_input", {
  tst <- ghcm_test(resid_X_on_Z_multivariate, resid_Y_on_Z_irregular,  Y_limits=c(0,2))
  expect_is(tst, "ghcm")
})

test_that("ghcm_runs_with_two_regular_inputs", {
  tst <- ghcm_test(resid_X_on_Z_multivariate, resid_Y_on_Z_scalar)
  expect_is(tst, "ghcm")
})

test_that("ghcm_fails_when_sample_sizes_differ_irregular", {
  resid_X_on_Z_irregular_half <- subset(resid_X_on_Z_irregular, .obs <= n/2)

  expect_error({
    ghcm_test(resid_X_on_Z_irregular_half, resid_Y_on_Z_irregular)
  })
})

test_that("ghcm_fails_when_sample_sizes_differ_irregular_and_regular", {
  resid_X_on_Z_multivariate_half <- resid_X_on_Z_multivariate[1:(n/2),]

  expect_error({
    ghcm_test(resid_X_on_Z_multivariate_half, resid_Y_on_Z_irregular)
  })
})

test_that("ghcm_fails_when_bad_X_limits_single_number", {

  expect_error({
    ghcm_test(resid_X_on_Z_irregular, resid_Y_on_Z_irregular, 0, c(0,2))
  })
})

test_that("ghcm_fails_when_bad_X_limits_bad_interval", {

  expect_error({
    ghcm_test(resid_X_on_Z_irregular, resid_Y_on_Z_irregular, c(1,0), c(0,2))
  })
})

test_that("ghcm_fails_when_bad_X_limits_single_number", {

  expect_error({
    ghcm_test(resid_X_on_Z_irregular, resid_Y_on_Z_irregular, c(0,1), 0)
  })
})

test_that("ghcm_fails_when_bad_X_limits_bad_interval", {

  expect_error({
    ghcm_test(resid_X_on_Z_irregular, resid_Y_on_Z_irregular, c(0,2), c(1,0))
  })
})

test_that("ghcm_fails_when_bad_resid_X", {
  bad_resid_X_on_Z_irregular <- list(.obs=1, .index=0.5, .value=1)

  expect_error({
    ghcm_test(bad_resid_X_on_Z_irregular, resid_Y_on_Z_irregular, c(0,2), c(0,1))
  })
})

test_that("ghcm_fails_when_bad_resid_X_wrong_columns", {
  bad_resid_X_on_Z_irregular <- data.frame(x=1, y=0.5, z=1)

  expect_error({
    ghcm_test(bad_resid_X_on_Z_irregular, resid_Y_on_Z_irregular, c(0,2), c(0,1))
  })
})


test_that("ghcm_fails_when_bad_resid_Y", {
  bad_resid_Y_on_Z_irregular <- list(.obs=1, .index=0.5, .value=1)

  expect_error({
    ghcm_test(resid_X_on_Z_irregular, bad_resid_Y_on_Z_irregular, c(0,2), c(0,1))
  })
})

test_that("ghcm_fails_when_bad_resid_Y_wrong_columns", {
  bad_resid_Y_on_Z_irregular <- data.frame(x=1, y=0.5, z=1)

  expect_error({
    ghcm_test(resid_X_on_Z_irregular, bad_resid_Y_on_Z_irregular, c(0,2), c(0,1))
  })
})

## Test gibberish
## Test data frame with wrong columns


