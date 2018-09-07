#' observed Pharmacometrics PK/PD Data
#'
#' A dataset containing concentration-time data with the
#' given dose and some subject characteristics to help in the app exploration.
#'
#' @format A data frame with 600 rows and 10 variables
#' \describe{
#'   \item{ID}{Subject Identifier, an integer from 57 to 186}
#'   \item{Time}{Time of sample measurement, in hours}
#'   \item{DV}{Dependent variable value}
#'   \item{DVTYPE}{type of Dependent variable with values "CObs","CMObs" or "EObs"}
#'   \item{REPL}{Simulation Repliate number, an integer from 0 to 50}
#' }
#' @source ASCPT Challenge competition
"simdata"