get_example_dir <- function() {
  example_dir <- system.file("extdata", package = "CellAgeR")
  if (!nzchar(example_dir)) {
    example_dir <- file.path("inst", "extdata")
  }
  example_dir
}

test_that("Walker example data produces predictions for packaged clocks", {
  clock_coefficients <- getClockCoefficients()

  expect_true(all(c("NeuronClock", "OligoClock") %in% names(clock_coefficients)))

  load(file.path(get_example_dir(), "Walker.Rdata"))

  clocks <- list(
    Neuron = clock_coefficients$NeuronClock,
    Oligo = clock_coefficients$OligoClock
  )

  cell_types <- c("NeuNPos", "Sox10Pos", "Total")

  for (clock_name in names(clocks)) {
    for (cell_type in cell_types) {
      idx <- which(phenoDf$CellType == cell_type)
      p_sub <- phenoDf[idx, ]
      b_sub <- as.matrix(betaMat[, idx, drop = FALSE])

      res <- predictCTSAge(b_sub, clocks[[clock_name]], age = p_sub$Age)

      expect_s3_class(res$predictions, "data.frame")
      expect_equal(nrow(res$predictions), ncol(b_sub))
      expect_true(all(c("SampleID", "DNAmAge") %in% colnames(res$predictions)))
    }
  }
})