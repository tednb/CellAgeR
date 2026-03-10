#' @title Access Packaged Clock Coefficients
#' @description Returns the packaged list of pre-trained cell-type specific
#'     clock coefficient tables included with CellAgeR.
#'
#' @return A named list of data frames. Each data frame contains the
#'     coefficients for one clock model.
#' @export
getClockCoefficients <- function() {
  clock_coefficients
}