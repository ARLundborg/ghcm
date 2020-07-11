ghcm_class_constructor <- function(test_statistic, p, limiting_covariance,
                       asymptotic_distribution, alpha) {
  #' Class constructor for the \code{ghcm} class.
  #'
  #' @param test_statistic Desc
  #' @param p Desc
  #' @param limiting_covariance Desc
  #' @param asymptotic_distribution Desc
  #' @param alpha Desc
  #'
  #' @export
  structure(class = "ghcm", list(
    test_statistic = test_statistic,
    p = p,
    limiting_covariance = limiting_covariance,
    asymptotic_distribution = asymptotic_distribution,
    alpha = alpha,
    reject = (p <= alpha))
    )
}

ghcm_test <- function(resid_X_on_Z, resid_Y_on_Z, b=10000,
                      X_is_multivariate=FALSE, Y_is_multivariate=FALSE,
                      fpca_method=NULL, alpha=0.05, ...) {
  #' Perform the GHCM using the residuals from regressing X on Z (resid_f)
  #'  and regressing Y on Z (resid_g).
  #'
  #'
  #' Desc
  #'
  #' @param resid_X_on_Z,resid_Y_on_Z This is a description. Bladebla
  #' @param X_is_multivariate,Y_is_multivariate This is a description. Bladebla
  #' @param b This is a description
  #' @param fpca_method This is a description. Bladebla
  #' @param ... Additional arguments to be passed to the fpca_method.
  #' @param alpha This is a description
  #'
  #' @export

  resid_X_on_Z <- as.matrix(resid_X_on_Z)
  resid_Y_on_Z <- as.matrix(resid_Y_on_Z)

  # Add n check

  # Add dim_X, dim_Y warning when not multivariate (i.e.
  # if d==1, always mult. if, say, d<10, print warning)

  if (!X_is_multivariate) {
    resid_X_on_Z <- fpca_method(resid_X_on_Z, ...)
  }
  if (!Y_is_multivariate) {
    resid_X_on_Z <- fpca_method(resid_X_on_Z, ...)
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
  asymptotic_dist <- rowSums(MASS::mvrnorm(b, rep(0, dim_X * dim_Y), cov_est)^2)
  p_value <- mean(asymptotic_dist > test_statistic)

  ghcm_class_constructor(test_statistic, p_value, cov_est, asymptotic_dist, alpha)
}



plot.ghcm <- function(x, ...) {
  #'
  #'
  #' @export
  if (x$p < 0.05) {
    xlim <- c(0, x$estimate * 1.25)
  }
  else {
    xlim <- c(0, max(x$asymptotic_dist))
  }
  graphics::hist(x$asymptotic_dist, prob = TRUE, xlim = xlim,
                 xlab = "Test statistic",
                 main = "Histogram of asymptotic test distribution", ...)
  graphics::abline(v = x$test_statistic, col = 2, lty = 2, lwd = 2, )
}

print.ghcm <- function(x, digits=getOption("digits"), ...) {
  #'
  #'
  #' @export
  cat("H0: X _||_ Y | Z, p:", format(x$p, digits = digits))
  cat("\n")
  invisible(x)
}
