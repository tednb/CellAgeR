#' @title Select Cell-Type Specific Features via SWRF
#' @description Identifies features associated with a specific phenotype 
#'     (continuous or categorical) in a cell-type-specific manner using a 
#'     Sliding Window Random Forest (SWRF) approach.
#'
#' @param targetCellType A string specifying the target cell type column name
#'     in `cellFractions`.
#' @param candidateFeatures A character vector of pre-selected feature names
#'     (e.g., DEGs, DMCs) to evaluate.
#' @param dataMatrix A numeric matrix of omics data (features in rows,
#'     samples in columns).
#' @param cellFractions A data frame of cell type fractions (samples in rows,
#'     cell types in columns).
#' @param responseVector A vector representing the target phenotype. Can be
#'     numeric (for regression) or factor/character (for classification).
#' @param seed Integer for random seed. Defaults to 42.
#' @param windowSize Number of samples per sliding window. Defaults to 100.
#' @param stepSize Number of samples to shift the window. Defaults to 10.
#' @param nTrees Number of trees for the Random Forest model. Defaults to 1000.
#' @param numThreads Number of CPU threads for ranger. Defaults to 20.
#'
#' @return A data frame containing Spearman correlation results for each
#'     feature.
#'
#' @examples
#' set.seed(1)
#' data_mat <- matrix(
#'   rnorm(20 * 24),
#'   nrow = 20,
#'   dimnames = list(paste0("gene", 1:20), paste0("Sample", 1:24))
#' )
#' cell_fracs <- data.frame(
#'   Neuron = seq(0.1, 0.9, length.out = 24),
#'   Astro = seq(0.9, 0.1, length.out = 24)
#' )
#' response <- 50 + 3 * cell_fracs$Neuron + rnorm(24, sd = 0.2)
#'
#' swrf_res <- SWRF(
#'   targetCellType = "Neuron",
#'   candidateFeatures = rownames(data_mat)[1:10],
#'   dataMatrix = data_mat,
#'   cellFractions = cell_fracs,
#'   responseVector = response,
#'   seed = 1,
#'   windowSize = 8,
#'   stepSize = 4,
#'   nTrees = 50,
#'   numThreads = 1
#' )
#'
#' head(swrf_res)
#' @export
SWRF <- function(targetCellType, candidateFeatures, dataMatrix,
                 cellFractions, responseVector, seed = 42,
                 windowSize = 100, stepSize = 10,
                 nTrees = 1000, numThreads = 20) {
  if (!is.matrix(dataMatrix) || !is.numeric(dataMatrix)) {
    stop("dataMatrix must be a numeric matrix.")
  }
  if (is.null(rownames(dataMatrix))) {
    stop("dataMatrix must have row names corresponding to feature IDs.")
  }
  if (!is.data.frame(cellFractions)) {
    stop("cellFractions must be a data frame.")
  }
  if (nrow(cellFractions) != ncol(dataMatrix)) {
    stop("nrow(cellFractions) must match ncol(dataMatrix).")
  }
  if (length(responseVector) != ncol(dataMatrix)) {
    stop("Length of responseVector must match ncol(dataMatrix).")
  }
  if (length(targetCellType) != 1L || !targetCellType %in% colnames(cellFractions)) {
    stop("targetCellType must name one column present in cellFractions.")
  }
  if (!is.numeric(windowSize) || length(windowSize) != 1L || windowSize < 2) {
    stop("windowSize must be a single numeric value greater than or equal to 2.")
  }
  if (!is.numeric(stepSize) || length(stepSize) != 1L || stepSize < 1) {
    stop("stepSize must be a single numeric value greater than or equal to 1.")
  }
  if (nrow(cellFractions) < windowSize) {
    stop("Total number of samples is less than windowSize.")
  }

  targetFracs <- cellFractions[[targetCellType]]
  sortedIdx <- order(targetFracs)

  sortedData <- dataMatrix[, sortedIdx, drop = FALSE]
  sortedResponse <- responseVector[sortedIdx]
  sortedFracs <- targetFracs[sortedIdx]

  if (is.character(sortedResponse)) {
    sortedResponse <- as.factor(sortedResponse)
  }
  isClassification <- is.factor(sortedResponse)

  featuresForRf <- intersect(rownames(sortedData), candidateFeatures)
  if (length(featuresForRf) == 0) {
    stop("No overlapping features between dataMatrix and candidateFeatures.")
  }

  dataSubset <- sortedData[featuresForRf, , drop = FALSE]
  windowStarts <- seq(1, ncol(dataSubset) - windowSize + 1, by = stepSize)
  numWindows <- length(windowStarts)

  if (isClassification) {
    message(sprintf("Running SWRF (Classification) across %d windows...", numWindows))
  } else {
    message(sprintf("Running SWRF (Regression) across %d windows...", numWindows))
  }

  impMatrix <- matrix(NA_real_, nrow = length(featuresForRf), ncol = numWindows)
  rownames(impMatrix) <- featuresForRf
  meanFracVec <- numeric(numWindows)

  for (i in seq_len(numWindows)) {
    startIdx <- windowStarts[i]
    windowIdx <- startIdx:(startIdx + windowSize - 1)

    rfModel <- ranger::ranger(
      x = t(dataSubset[, windowIdx, drop = FALSE]),
      y = sortedResponse[windowIdx],
      num.trees = nTrees,
      importance = "permutation",
      num.threads = numThreads,
      seed = seed
    )

    impMatrix[, i] <- rfModel$variable.importance[featuresForRf]
    meanFracVec[i] <- mean(sortedFracs[windowIdx], na.rm = TRUE)
  }

  message("Calculating Spearman correlations...")
  corEst <- numeric(length(featuresForRf))
  pVals <- numeric(length(featuresForRf))

  for (j in seq_along(featuresForRf)) {
    impScores <- impMatrix[j, ]

    if (stats::sd(impScores, na.rm = TRUE) > 0) {
      testRes <- stats::cor.test(
        impScores,
        meanFracVec,
        method = "spearman",
        exact = FALSE
      )
      corEst[j] <- unname(testRes$estimate)
      pVals[j] <- testRes$p.value
    } else {
      corEst[j] <- 0
      pVals[j] <- 1
    }
  }

  resultDf <- data.frame(
    Feature = featuresForRf,
    Correlation = corEst,
    P_Value = pVals,
    stringsAsFactors = FALSE
  )
  resultDf$P_Adj <- stats::p.adjust(resultDf$P_Value, method = "BH")

  resultDf <- resultDf[order(resultDf$P_Adj, -abs(resultDf$Correlation)), ]
  rownames(resultDf) <- NULL

  message("SWRF feature selection complete.")
  resultDf
}