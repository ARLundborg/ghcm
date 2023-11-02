ghcm_class_constructor <- function(test_statistic, p, alpha) {
  #' Class constructor for the \code{ghcm} class.
  #'
  #' @param test_statistic Positive numeric. The value of the GHCM test
  #' statistic.
  #' @param p Numeric in the unit interval. The estimated p-value of the
  #'  statistic.
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
    alpha = alpha,
    reject = (p <= alpha))
    )
}


ghcm_test <- function(resid_X_on_Z, resid_Y_on_Z, X_limits=NULL, Y_limits=NULL, alpha=0.05) {
  #' Conditional Independence Test using the GHCM
  #'
  #'
  #' Test whether X is independent of Y given Z using the Generalised Hilbertian
  #'  Covariance Measure. The function is applied to residuals from regressing
  #'  each of X and Y on Z respectively. Its validity is contingent on the performance
  #'  of the regression methods. For a more in-depth explanation see the package
  #'  vignette or the paper mentioned in the references.
  #'
  #'
  #' @param resid_X_on_Z,resid_Y_on_Z Residuals from regressing X (Y) on Z with
  #'  a suitable regression method. If X (Y) is uni- or multivariate or
  #'  functional on a constant, fixed grid, the residuals should be supplied as
  #'  a vector or matrix with no missing values. If instead X (Y) is functional
  #'  and observed on varying grids or with missing values, the residuals should
  #'  be supplied as a "melted" data frame with
  #'   \describe{
  #'       \item{.obs}{Integer indicating which curve the row corresponds to.}
  #'       \item{.index}{Function argument that the curve is evaluated at.}
  #'       \item{.value}{Value of the function.}
  #'   }
  #'   Note that in the irregular case, a minimum of 4 observations per curve is
  #'   required.
  #' @param X_limits,Y_limits The minimum and maximum values of the function
  #' argument of the X (Y) curves. Ignored if X (Y) is not functional.
  #' @param alpha Numeric in the unit interval. Significance level of the test.
  #'
  #' @return An object of class \code{ghcm} containing:
  #'   \describe{
  #'     \item{\code{test_statistic}}{Numeric, test statistic of the test.}
  #'     \item{\code{p}}{Numeric in the unit interval, estimated p-value of
  #'      the test.}
  #'     \item{\code{alpha}}{Numeric in the unit interval, significance level
  #'      of the test.}
  #'     \item{\code{reject}}{\code{TRUE} if \code{p} < \code{alpha}, \code{FALSE} otherwise.}
  #'   }
  #'
  #' @references
  #'  Please cite the following paper: Anton Rask Lundborg, Rajen D. Shah and
  #'  Jonas Peters: "Conditional Independence Testing in Hilbert Spaces with
  #'  Applications to Functional Data Analysis" Journal of the Royal Statistical
  #'  Society: Series B (Statistical Methodology) 2022 \doi{10.1111/rssb.12544}.
  #'
  #' @examples
  #' if (require(refund)) {
  #'   set.seed(1)
  #'   data(ghcm_sim_data)
  #'   grid <- seq(0, 1, length.out = 101)
  #'
  #' # Test independence of two scalars given a functional variable
  #'
  #'   m_1 <- pfr(Y_1 ~ lf(Z), data=ghcm_sim_data)
  #'   m_2 <- pfr(Y_2 ~ lf(Z), data=ghcm_sim_data)
  #'   ghcm_test(resid(m_1), resid(m_2))
  #'
  #' # Test independence of a regularly observed functional variable and a
  #' # scalar variable given a functional variable
  #'   \donttest{
  #'     m_X <- pffr(X ~ ff(Z), data=ghcm_sim_data, chunk.size=31000)
  #'     ghcm_test(resid(m_X), resid(m_1))
  #'   }
  #' # Test independence of two regularly observed functional variables given
  #' # a functional variable
  #'   \donttest{
  #'      m_W <- pffr(W ~ ff(Z), data=ghcm_sim_data, chunk.size=31000)
  #'     ghcm_test(resid(m_X), resid(m_W))
  #'   }
  #'
  #'
  #'   data(ghcm_sim_data_irregular)
  #'   n <- length(ghcm_sim_data_irregular$Y_1)
  #'   Z_df <- data.frame(.obs=1:n)
  #'   Z_df$Z <- ghcm_sim_data_irregular$Z
  #' # Test independence of an irregularly observed functional variable and a
  #' # scalar variable given a functional variable
  #'   \donttest{
  #'     m_1 <- pfr(Y_1 ~ lf(Z), data=ghcm_sim_data_irregular)
  #'     m_X <- pffr(X ~ ff(Z), ydata = ghcm_sim_data_irregular$X,
  #'     data=Z_df, chunk.size=31000)
  #'     ghcm_test(resid(m_X), resid(m_1), X_limits=c(0, 1))
  #'  }
  #' # Test independence of two irregularly observed functional variables given
  #' # a functional variable
  #'   \donttest{
  #'     m_W <- pffr(W ~ ff(Z), ydata = ghcm_sim_data_irregular$W,
  #'     data=Z_df, chunk.size=31000)
  #'     ghcm_test(resid(m_X), resid(m_W), X_limits=c(0, 1), Y_limits=c(0, 1))
  #'  }
  #' }
  #'
  #' @export


  if(is.data.frame(resid_X_on_Z)) {

    if (!is.numeric(X_limits)) {
      stop("X_limits needs to be an interval when X is functional.")
    } else if(length(X_limits) != 2) {
      stop("X_limits needs to be an interval when X is functional.")
    } else if(X_limits[1] > X_limits[2]) {
      stop("X_limits needs to be an interval when X is functional.")
    }

    if(!all(c(".obs", ".value", ".index") %in% colnames(resid_X_on_Z))) {
      stop("resid_Y_on_Z does not contain the required columns '.obs', '.index' and '.value'.")
    }

    resid_X_on_Z_split <- split(resid_X_on_Z, resid_X_on_Z$.obs)
    n <- length(resid_X_on_Z_split)

    resid_X_on_Z_splines <- lapply(resid_X_on_Z_split, function(x)
                              splines::interpSpline(x$.index, x$.value))
    resid_X_on_Z_prods <- inner_product_matrix_splines(resid_X_on_Z_splines, X_limits[1], X_limits[2])

  } else if(is.numeric(resid_X_on_Z)) {
    resid_X_on_Z <- as.matrix(resid_X_on_Z)
    n <- dim(resid_X_on_Z)[1]

    centered_resid_X_on_Z <- sweep(resid_X_on_Z, 2, colMeans(resid_X_on_Z))
    resid_X_on_Z_prods <- tcrossprod(centered_resid_X_on_Z)
  } else{
    stop("resid_X_on_Z is an unrecognised data type.")
  }

  if(is.data.frame(resid_Y_on_Z)) {

    if (!is.numeric(Y_limits)) {
      stop("Y_limits needs to be an interval when X is functional.")
    } else if(length(Y_limits) != 2) {
      stop("Y_limits needs to be an interval when X is functional.")
    } else if(Y_limits[1] > Y_limits[2]) {
      stop("Y_limits needs to be an interval when X is functional.")
    }

    if(!all(c(".obs", ".value", ".index") %in% colnames(resid_Y_on_Z))) {
      stop("resid_Y_on_Z does not contain the required columns '.obs', '.index' and '.value'.")
    }

    resid_Y_on_Z_split <- split(resid_Y_on_Z, resid_Y_on_Z$.obs)
    if(length(resid_Y_on_Z_split) != n) {
      stop("The sample sizes of the X residuals and Y residuals differ.")
    }

    resid_Y_on_Z_splines <- lapply(resid_Y_on_Z_split, function(x)
      splines::interpSpline(x$.index, x$.value))
    resid_Y_on_Z_prods <- inner_product_matrix_splines(resid_Y_on_Z_splines, Y_limits[1], Y_limits[2])

  } else if(is.numeric(resid_Y_on_Z)) {
    resid_Y_on_Z <- as.matrix(resid_Y_on_Z)
    if(dim(resid_Y_on_Z)[1] != n) {
      stop("The sample sizes of the X residuals and Y residuals differ.")
    }

    centered_resid_Y_on_Z <- sweep(resid_Y_on_Z, 2, colMeans(resid_Y_on_Z))
    resid_Y_on_Z_prods <- tcrossprod(centered_resid_Y_on_Z)
  } else{
    stop("resid_Y_on_Z is an unrecognised data type.")
  }

  inner_product_mat <- resid_X_on_Z_prods*resid_Y_on_Z_prods

  n <- dim(inner_product_mat)[1]

  test_statistic <- 1/n * sum(inner_product_mat)
  scaling_mat <- matrix(1/n, n, n)
  cov_mat <- 1/(n - 1) * (inner_product_mat - scaling_mat %*%
                            inner_product_mat - inner_product_mat %*% scaling_mat +
                            scaling_mat %*% inner_product_mat %*% scaling_mat)
  eig_vals <- pmax(eigen(cov_mat, only.values = TRUE, symmetric = TRUE)$values[-n],
                   0)
  p <- tryCatch(CompQuadForm::imhof(test_statistic, eig_vals,
                                    epsrel = .Machine$double.eps, epsabs = .Machine$double.eps)$Qq,
                warning = function(w) return(0))
  ghcm_class_constructor(test_statistic, p, alpha)
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
