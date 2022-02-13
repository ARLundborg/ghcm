#' GHCM simulated data
#'
#' A simulated dataset containing a combination of functional and scalar
#'  variables. Y_1 and Y_2 are scalar random variables and are both functions
#'  of Z. X, Z and W are functional, Z is a function of X and W is a function of Z.
#'
#' In \code{ghcm_sim_data} the functional variables each consists of 101 observations on
#'  an equidistant grid on [0, 1].
#'
#' In \code{ghcm_sim_data_irregular} the functional variables X and W are instead only observed on
#' a subsample of the original equidistant grid.
#'
#' @format \code{ghcm_sim_data} is a data frame with 500 rows of 5 variables:
#' \describe{
#'   \item{Y_1}{Numeric vector.}
#'   \item{Y_2}{Numeric vector.}
#'   \item{Z}{500 x 101 matrix.}
#'   \item{X}{500 x 101 matrix.}
#'   \item{W}{500 x 101 matrix.}
#' }
#'
#' @format \code{ghcm_sim_data_irregular} is a list with 5 elements:
#' \describe{
#'   \item{Y_1}{Numeric vector.}
#'   \item{Y_2}{Numeric vector.}
#'   \item{Z}{500 x 101 matrix.}
#'   \item{X}{A data frame with
#'     \describe{
#'       \item{.obs}{Integer between 1 and 500 indicating which curve the row corresponds to.}
#'       \item{.index}{Function argument that the curve is evaluated at.}
#'       \item{.value}{Value of the function.}
#'     }}
#'   \item{W}{A data frame with
#'     \describe{
#'       \item{.obs}{Integer between 1 and 500 indicating which curve the row corresponds to.}
#'       \item{.index}{Function argument that the curve is evaluated at.}
#'       \item{.value}{Value of the function.}
#'     }}
#' }
#'
#' @source The generation script can be found in the \code{data-raw} folder of
#' the package.
#' @rdname ghcm_sim_data
"ghcm_sim_data"

#' @name ghcm_sim_data
#' @format NULL
#' @rdname ghcm_sim_data
"ghcm_sim_data_irregular"

