ghcm_class_constructor <- function(test_statistic, p, dim, cov, samples,
                                   alpha) {
  #' Class constructor for the \code{ghcm} class.
  #'
  #' @param test_statistic Positive numeric. The value of the GHCM test
  #' statistic.
  #' @param p Numeric in the unit interval. The estimated p-value of the
  #'  statistic.
  #' @param dim Positive integer. The dimension of the truncated Gaussian limit
  #'  distribution.
  #' @param cov \code{dim} x \code{dim} covariance matrix. The covariance of the
  #'  truncated Gaussian limit distribution.
  #' @param samples Numeric vector. Samples from the estimated limiting
  #'  distribution.
  #' @param alpha Numeric in the unit interval. The significance level of the
  #'  test.
  #' @keywords internal
  #' @return An object of class \code{ghcm} with the properties listed above.
  #'
  #' @noRd
  #'
  structure(class = "ghcm", list(
    test_statistic = test_statistic,
    p = p,
    dim = as.integer(dim),
    cov = cov,
    samples = samples,
    alpha = alpha,
    reject = (p <= alpha))
    )
}

ghcm_test <- function(resid_X_on_Z, resid_Y_on_Z,
               X_grid=NULL,
               Y_grid=NULL,
               fpca_method="fpca.sc", b=10000, alpha=0.05, ...) {
  #' Conditional Independence Test using the GHCM
  #'
  #'
  #' Testing X independent of Y given Z using the Generalised Hilbertian
  #'  Covariance Measure. The function is applied to residuals from regressing X on Z
  #'  and regressing Y on Z and its validity is contingent on the performance
  #'  of the regression methods.
  #'
  #' @param resid_X_on_Z,resid_Y_on_Z Numeric vectors or matrices. Residuals
  #'  when regressing X (Y) on Z with a suitable regression method.
  #' @param X_grid,Y_grid Numeric vectors or NA. The grid of values that X (Y)
  #'  is observed on. When set to NULL (the default) an equidistant grid on the
  #'  unit interval is created. If NA, X (Y) is assumed to not be a functional
  #'  random variable.
  #' @param fpca_method String or function. If a string is given, will search
  #'  the refund package for a function with the given name. If a function is
  #'  given it must take a data matrix and a grid as input and return a matrix
  #'  with the same number of rows as the input and the coordinates of the input
  #'  in its FPCA basis as each row. Extra arguments to the fpca function
  #'  are supplied with \code{...}. Currently supported refund fpca functions
  #'  are \code{fpca.sc} (the default), \code{fpca.ssvd} and \code{fpca.face}.
  #' @param b Positive integer. The number of Monte Carlo realisations from the
  #'  estimated limiting distribution to estimate the p-value.
  #' @param alpha Numeric in the unit interval. Significance level of the test.
  #' @param ... Additional arguments to be passed to the fpca_method.
  #'
  #' @return An object of class \code{ghcm} containing:
  #'   \describe{
  #'     \item{\code{test_statistic}}{Numeric, test statistic of the test.}
  #'     \item{\code{p}}{Numeric in the unit interval, estimated p-value of
  #'      the test.}
  #'     \item{\code{dim}}{Positive integer, the dimension of the truncated
  #'      limiting Gaussian.}
  #'     \item{\code{cov}}{\code{dim} x \code{dim} matrix, estimated covariance
  #'      of the truncated limiting Gaussian.}
  #'     \item{\code{samples}}{Numeric vector, samples of the Hilbert-Schmidt
  #'      norm of the estimated truncated limiting Gaussian.}
  #'     \item{\code{alpha}}{Numeric in the unit interval, significance level
  #'      of the test.}
  #'   }
  #'
  #' @references
  #'  Please cite the following paper: Anton Rask Lundborg, Rajen D. Shah and
  #'  Jonas Peters: "Conditional Independence Testing in Hilbert Spaces with
  #'  Applications to Functional Data Analysis" https://arxiv.org/abs/2101.07108
  #'
  #' @examples
  #' library(refund)
  #' set.seed(1)
  #' data(ghcm_sim_data)
  #' grid <- seq(0, 1, length.out = 101)
  #'
  #' # Test independence of two scalars given a functional variable
  #'
  #' m_1 <- pfr(Y_1 ~ lf(Z), data=ghcm_sim_data)
  #' m_2 <- pfr(Y_2 ~ lf(Z), data=ghcm_sim_data)
  #' ghcm_test(resid(m_1), resid(m_2), X_grid = NA, Y_grid = NA )
  #'
  #' # Test independence of a functional variable and a scalar variable given a
  #' # functional variable
  #' \donttest{
  #' m_X <- pffr(X ~ ff(Z), data=ghcm_sim_data, chunk.size=31000)
  #' ghcm_test(resid(m_X), resid(m_1), X_grid = grid, Y_grid = NA )
  #'}
  #' # Test independence of two functional variables given a functional variable
  #' \donttest{
  #' m_W <- pffr(W ~ ff(Z), data=ghcm_sim_data, chunk.size=31000)
  #' ghcm_test(resid(m_X), resid(m_W), X_grid = grid, Y_grid = grid)
  #'}
  #' @export


  resid_X_on_Z <- as.matrix(resid_X_on_Z)
  resid_Y_on_Z <- as.matrix(resid_Y_on_Z)

  if (dim(resid_X_on_Z)[1] != dim(resid_Y_on_Z)[1]) {
    stop("Error: The sample sizes of the X residuals and Y residuals differ.")
  }

  if (is.null(X_grid)) {
    X_grid <- seq(0, 1, length.out = dim(resid_X_on_Z)[2])
    X_manually_set <- FALSE
  }
  else {
    X_manually_set <- TRUE
  }

  if (is.null(Y_grid)) {
    Y_grid <- seq(0, 1, length.out = dim(resid_Y_on_Z)[2])
    Y_manually_set <- FALSE
  }
  else {
    Y_manually_set <- TRUE
  }


  if (is.character(fpca_method)) {
    tryCatch(
      expr = refund_fpca_function <- utils::getFromNamespace(fpca_method,
                                                             "refund"),
      error = function(e) {
        stop(paste0("Error: Could not find '", fpca_method,
                    "' in the refund package."))
      })
    fpca_method <- function(x, grid, ...) {
      refund_fpca_function(x, argvals = grid, ...)$scores
      }
  }
  else if (is.function(fpca_method)) {
    fpca_method <- fpca_method
  }
  else {
    stop("Error: fpca_method should be a function or a character string.")
  }

  X_not_functional <- all(is.na(X_grid))

  if (X_not_functional | length(X_grid) == 1) {
    # Do not do FPCA
  }
  else {
    if (!X_manually_set) {
      warning("Warning: X is assumed to have been observed on an equidistant
      grid. If this is not the case or if you want to suppress the warning set
              X_grid.")
    }
    resid_X_on_Z <- fpca_method(resid_X_on_Z, X_grid, ...)
  }

  Y_not_functional <- all(is.na(Y_grid))

  if (Y_not_functional | length(Y_grid) == 1) {
    # Do not do FPCA
  }
  else {
    if (!Y_manually_set) {
      warning("Warning: Y is assumed to have been observed on an equidistant
      grid. If this is not the case or if you want to suppress the warning set
              Y_grid.")
    }
    resid_Y_on_Z <- fpca_method(resid_Y_on_Z, Y_grid, ...)
  }


  n <- dim(resid_X_on_Z)[1]
  dim_X <- dim(resid_X_on_Z)[2]
  dim_Y <- dim(resid_Y_on_Z)[2]

  outer_products <- t(apply(cbind(resid_X_on_Z, resid_Y_on_Z), 1,
                            function(x) {
                              x[1:dim_X] %o% x[(dim_X + 1):(dim_X + dim_Y)]
                            }
  ))

  if (dim(outer_products)[1] == 1) {
    outer_products <- t(outer_products)
  }

  test_statistic <- sum((sqrt(n) * colMeans(outer_products))^2)
  cov_est <- stats::cov(outer_products)
  limit_dim <- dim_X * dim_Y
  samples <- rowSums(MASS::mvrnorm(b, rep(0, limit_dim), cov_est)^2)
  p <- (1+sum(samples > test_statistic))/(1+b)

  ghcm_class_constructor(test_statistic, p, limit_dim, cov_est, samples, alpha)
}

plot.ghcm <- function(x, bw="SJ", ...) {
  #' Plotting function for the \code{ghcm}-class.
  #'
  #' Plots the observed test statistic of a performed GHCM test together with a
  #' density estimate of the estimated asymptotic distribution of the test
  #' statistic under the null.
  #'
  #' @param x \code{ghcm}-object. The result of running a \code{ghcm_test}.
  #' @param bw a string or numeric. If \code{bw} is a string, it should specify
  #' a bandwidth method for the \code{density} function in the \pkg{stats}
  #' package. If \code{bw} is a numeric, it will be used as the bandwidth in
  #' the call to \code{density}.
  #' @param ... additional plotting parameters.
  #'
  #' @return
  #' None.
  #'
  #' @examples
  #' set.seed(1)
  #' library(refund)
  #' data(ghcm_sim_data)
  #'
  #' # Test independence of two scalars given a functional variable
  #'
  #' m_1 <- pfr(Y_1 ~ lf(Z), data=ghcm_sim_data)
  #' m_2 <- pfr(Y_2 ~ lf(Z), data=ghcm_sim_data)
  #' test <- ghcm_test(resid(m_1), resid(m_2), X_grid = NA, Y_grid = NA )
  #' plot(test)
  #'
  #' @export
  #'

  density_estimate <- stats::density(x$samples, from=0, bw=bw)
  if(max(x$samples) >= x$test_statistic) {
    graphics::plot(density_estimate,
                   xlab = "Test statistic",
                   main = "Density estimate of asymptotic test distribution", ...)
  }
  else {
    graphics::plot(density_estimate,
                   xlab = "Test statistic",
                   main = "Density estimate of asymptotic test distribution", xlim=c(0, x$test_statistic))
  }
  graphics::rug(x$samples)
  graphics::abline(v = x$test_statistic, col = 2, lty = 2, lwd = 2)
}

print.ghcm <- function(x, digits=getOption("digits"), ...) {
  #' @export
  cat("H0: X _||_ Y | Z, p:", format(x$p, digits = digits))
  cat("\n")
  if(x$reject){
    cat("Rejected at",format(x$alpha*100, digits = digits),"% level")
    cat("\n")
  } else {
    cat("Not rejected at",format(x$alpha*100, digits = digits),"% level")
    cat("\n")
  }
  cat("Test statistic:", format(x$test_statistic, digits = digits))
  cat("\n")
  invisible(x)
}
