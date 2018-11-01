#' Simulated Pharmacometrics PK/PD Data
#'
#' A simulated dataset containing concentration-time data with the
#' type of dependent variable and replicate number.
#'
#' @format A data frame with 357700 rows and 5 variables
#' \describe{
#'   \item{ID}{Subject Identifier, an integer from 57 to 186}
#'   \item{TIME}{Time of sample measurement, in hours}
#'   \item{DV}{Dependent variable value}
#'   \item{DVTYPE}{type of Dependent variable with values "CObs","CMObs" or "EObs"}
#'   \item{REPL}{Simulation Replicate number, an integer from 0 to 50}
#' }
#' @source ASCPT Challenge competition
"simdata"