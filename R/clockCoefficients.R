#' @title Access Packaged Clock Coefficients
#' @description Returns the packaged list of pre-trained cell-type specific
#'     clock coefficient tables included with CellAgeR. Available clock names
#'     are conceptually structured by tissue and cell type:
#'     Brain (`NeuronClock`, `sortedNeuronClock`, `OligoClock`, `sortedOligoClock`);
#'     Epithelial (`BreastEpiClock`, `ColonEpiClock`, `KeraClock`, `LungEpiClock`,
#'     `ProstateEpiClock`, `salivaEpiClock`); Fibroblast (`BreastFibClock`,
#'     `ColonFibClock`, `LungFibClock`); Immune (`Cd4tClock`, `ColonImmunClock`,
#'     `LungImmunClock`, `MonoClock`, `salivaImmunClock`); and Liver (`HepClock`).
#'
#' @return A named list of data frames. Each data frame contains the
#'     coefficients for one clock model.
#' @examples
#' names(getClockCoefficients())
#' @export
getClockCoefficients <- function() {
  clock_coefficients
}