#' GHCM simulated data
#'
#' A simulated dataset containing a combination of functional and scalar
#'  variables. The functional variables each consists of 101 observations on
#'  an equidistant grid on [0, 1].
#'
#' Y_1 and Y_2 are scalar random variables and are both functions of Z. X, Z and
#'  W are functional, Z is a function of X and W is a function of Z.
#'
#'
#' @format A data frame with 500 rows of 5 variables:
#' \describe{
#'   \item{X}{500 x 101 matrix.}
#'   \item{Z}{500 x 101 matrix.}
#'   \item{W}{500 x 101 matrix. }
#'   \item{Y_1}{Numeric vector.}
#'   \item{Y_2}{Numeric vector.}
#' }
#' @source The generation script can be found in the \code{data-raw} folder of
#' the package.
#'
"ghcm_sim_data"
