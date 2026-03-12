test_that("calculateCSS returns expected columns", {
  set.seed(1)
  beta_mat <- matrix(
    runif(10 * 6, min = 0.05, max = 0.95),
    nrow = 10,
    dimnames = list(paste0("cg", 1:10), paste0("Sample", 1:6))
  )
  cell_types <- c("Neuron", "Neuron", "Neuron", "Astro", "Astro", "Astro")

  res <- calculateCSS(
    betaMatrix = beta_mat,
    cellTypes = cell_types,
    targetCellType = "Neuron"
  )

  expect_s3_class(res, "data.frame")
  expect_equal(
    colnames(res),
    c("Feature", "logFC", "P_Value", "FDR", "CSS")
  )
  expect_equal(nrow(res), nrow(beta_mat))
})

test_that("SWRF returns ranked feature statistics", {
  set.seed(1)
  data_mat <- matrix(
    rnorm(20 * 24),
    nrow = 20,
    dimnames = list(paste0("gene", 1:20), paste0("Sample", 1:24))
  )
  cell_fracs <- data.frame(
    Neuron = seq(0.1, 0.9, length.out = 24),
    Astro = seq(0.9, 0.1, length.out = 24)
  )
  response <- 50 + 3 * cell_fracs$Neuron + rnorm(24, sd = 0.2)

  res <- SWRF(
    targetCellType = "Neuron",
    candidateFeatures = rownames(data_mat)[1:10],
    dataMatrix = data_mat,
    cellFractions = cell_fracs,
    responseVector = response,
    seed = 1,
    windowSize = 8,
    stepSize = 4,
    nTrees = 50,
    numThreads = 1
  )

  expect_s3_class(res, "data.frame")
  expect_equal(
    colnames(res),
    c("Feature", "Correlation", "P_Value", "P_Adj")
  )
  expect_true(nrow(res) > 0)
})

test_that("trainWeightedClock returns a fitted clock", {
  set.seed(1)
  data_mat <- matrix(
    rnorm(30 * 20),
    nrow = 30,
    dimnames = list(paste0("cg", 1:30), paste0("Sample", 1:20))
  )
  response <- seq(30, 68, length.out = 20) + rnorm(20, sd = 1)
  candidate_features <- rownames(data_mat)[1:12]
  feature_weights <- setNames(seq(0.5, 1.5, length.out = 12), candidate_features)

  res <- trainWeightedClock(
    dataMatrix = data_mat,
    responseVector = response,
    candidateFeatures = candidate_features,
    featureWeights = feature_weights,
    alpha = 0.5,
    nfolds = 5,
    seed = 1
  )

  expect_type(res, "list")
  expect_true(all(c("cvModel", "bestLambda", "clockCoefficients") %in% names(res)))
  expect_s3_class(res$clockCoefficients, "data.frame")
  expect_true(all(c("Feature", "Coefficient") %in% colnames(res$clockCoefficients)))
})