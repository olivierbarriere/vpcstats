#' observed Pharmacometrics PK/PD Data
#'
#' A dataset containing concentration-time data with the
#' type of dependent variable.
#'
#' @format A data frame with 7154 rows and 4 variables
#' \describe{
#'   \item{ID}{Subject Identifier, an integer from 57 to 186}
#'   \item{TIME}{Time of sample measurement, in hours}
#'   \item{DV}{Dependent variable value}
#'   \item{DVTYPE}{type of Dependent variable with values "CObs","CMObs" or "EObs"}
#' }
#' @source ASCPT Challenge competition
"obsdata"