ghcm_class_constructor <- function(test_statistic, p, cov,
                                   alpha) {
  #' Class constructor for the \code{ghcm} class.
  #'
  #' @param test_statistic Positive numeric. The value of the GHCM test
  #' statistic.
  #' @param p Numeric in the unit interval. The estimated p-value of the
  #'  statistic.
  #' @param cov \code{dim} x \code{dim} covariance matrix. The covariance of the
  #'  truncated Gaussian limit distribution.
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
    cov = cov,
    alpha = alpha,
    reject = (p <= alpha))
    )
}

ghcm_test <- function(resid_X_on_Z, resid_Y_on_Z, alpha=0.05) {
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
  #' @param alpha Numeric in the unit interval. Significance level of the test.
  #'
  #' @return An object of class \code{ghcm} containing:
  #'   \describe{
  #'     \item{\code{test_statistic}}{Numeric, test statistic of the test.}
  #'     \item{\code{p}}{Numeric in the unit interval, estimated p-value of
  #'      the test.}
  #'     \item{\code{cov}}{matrix, estimated covariance
  #'      of the truncated limiting Gaussian.}
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
  #' ghcm_test(resid(m_1), resid(m_2))
  #'
  #' # Test independence of a functional variable and a scalar variable given a
  #' # functional variable
  #' \donttest{
  #' m_X <- pffr(X ~ ff(Z), data=ghcm_sim_data, chunk.size=31000)
  #' ghcm_test(resid(m_X), resid(m_1))
  #'}
  #' # Test independence of two functional variables given a functional variable
  #' \donttest{
  #' m_W <- pffr(W ~ ff(Z), data=ghcm_sim_data, chunk.size=31000)
  #' ghcm_test(resid(m_X), resid(m_W))
  #'}
  #' @export


  resid_X_on_Z <- as.matrix(resid_X_on_Z)
  resid_Y_on_Z <- as.matrix(resid_Y_on_Z)

  resid_X_on_Z <- sweep(resid_X_on_Z, 2, colMeans(resid_X_on_Z))
  resid_Y_on_Z <- sweep(resid_Y_on_Z, 2, colMeans(resid_Y_on_Z))

  if (dim(resid_X_on_Z)[1] != dim(resid_Y_on_Z)[1]) {
    stop("The sample sizes of the X residuals and Y residuals differ.")
  }

  n <- dim(resid_X_on_Z)[1]

  inner_product_mat <- tcrossprod(resid_X_on_Z)*tcrossprod(resid_Y_on_Z)

  test_statistic <- 1/n*sum(inner_product_mat)

  scaling_mat <- matrix(1/n, n, n)

  cov_mat <- 1/(n-1)*(inner_product_mat - scaling_mat %*% inner_product_mat - inner_product_mat %*% scaling_mat + scaling_mat %*% inner_product_mat %*% scaling_mat)

  eig_vals <- pmax(eigen(cov_mat, only.values = TRUE, symmetric=TRUE)$values[-n], 0)

  p <- tryCatch(
    CompQuadForm::imhof(test_statistic, eig_vals, epsrel = .Machine$double.eps,
                           epsabs = .Machine$double.eps)$Qq,
    warning=function(w) return(0))
  ghcm_class_constructor(test_statistic, p, cov_mat, alpha)
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
  invisible(x)
}
