get_example_dir <- function() {
  example_dir <- system.file("extdata", package = "CellAgeR")
  if (!nzchar(example_dir)) {
    example_dir <- file.path("inst", "extdata")
  }
  example_dir
}

test_that("HD striatum example supports age acceleration analysis", {
  clock_coefficients <- getClockCoefficients()
  load(file.path(get_example_dir(), "HDStriatum.Rdata"))

  res <- predictCTSAge(
    betaMatrix = as.matrix(betaMat),
    clock = clock_coefficients$NeuronClock,
    age = phenoDf$Age
  )

  expect_s3_class(res$predictions, "data.frame")
  expect_s3_class(res$ageAcceleration, "data.frame")
  expect_equal(nrow(res$ageAcceleration), ncol(betaMat))

  plot_df <- data.frame(
    Group = phenoDf$GroupRaw,
    AgeAccel = res$ageAcceleration$AgeAcceleration
  )

  plot_df <- na.omit(plot_df)
  plot_df <- plot_df[plot_df$AgeAccel >= -30, ]

  stat_res <- wilcox.test(AgeAccel ~ Group, data = plot_df, alternative = "greater")

  expect_type(stat_res$p.value, "double")
  expect_length(stat_res$p.value, 1)
  expect_true(is.finite(stat_res$p.value))
})