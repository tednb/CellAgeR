#' @title Access Packaged Clock Coefficients
#' @description Returns the packaged list of pre-trained cell-type specific
#'     clock coefficient tables included with CellAgeR. Available clock names
#'     are `NeuronClock`, `OligoClock`, `sortedNeuronClock`,
#'     `sortedOligoClock`, `Cd4tClock`, `MonoClock`, `ColonEpiClock`,
#'     `KeraClock`, `LungEpiClock`, and `HepClock`.
#'
#' @return A named list of data frames. Each data frame contains the
#'     coefficients for one clock model.
#' @examples
#' names(getClockCoefficients())
#' @export
getClockCoefficients <- function() {
  clock_coefficients
}