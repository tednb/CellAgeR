#' @title Train a Cell Type Specific Clock
#' @description Trains an Elastic Net model using pre-selected features and 
#'     custom feature weights (penalties). This function serves as the core 
#'     predictive model builder for generating custom or cell-type-specific clocks.
#'
#' @param dataMatrix A numeric matrix of DNAm data (features in rows,
#'     samples in columns).
#' @param responseVector A numeric vector representing the target phenotype 
#'     (e.g., chronological age). Length must match `ncol(dataMatrix)`.
#' @param candidateFeatures A character vector of features to include in the 
#'     model training. Must intersect with `rownames(dataMatrix)`.
#' @param featureWeights A named numeric vector of weights applied to the 
#'     `candidateFeatures`. Note: These are passed as `penalty.factor` to 
#'     glmnet, meaning a LOWER value indicates a MORE important feature 
#'     (less shrinkage). If NULL, defaults to unweighted (all 1s).
#' @param alpha The elastic net mixing parameter. Defaults to 0.5.
#' @param nfolds Number of folds for cross-validation. Defaults to 10.
#' @param seed Integer for random seed to ensure reproducible CV folds. 
#'     Defaults to 42.
#'
#' @return A list containing `cvModel`, `bestLambda`, and `clockCoefficients`.
#'   The `clockCoefficients` component is a data frame of non-zero coefficients,
#'   including the intercept when present.
#'
#' @examples
#' set.seed(1)
#' data_mat <- matrix(
#'   rnorm(30 * 20),
#'   nrow = 30,
#'   dimnames = list(paste0("cg", 1:30), paste0("Sample", 1:20))
#' )
#' response <- seq(30, 68, length.out = 20) + rnorm(20, sd = 1)
#' candidate_features <- rownames(data_mat)[1:12]
#' feature_weights <- setNames(seq(0.5, 1.5, length.out = 12), candidate_features)
#'
#' clock_res <- trainWeightedClock(
#'   dataMatrix = data_mat,
#'   responseVector = response,
#'   candidateFeatures = candidate_features,
#'   featureWeights = feature_weights,
#'   alpha = 0.5,
#'   nfolds = 5,
#'   seed = 1
#' )
#'
#' head(clock_res$clockCoefficients)
#' @export
trainWeightedClock <- function(dataMatrix, responseVector, candidateFeatures,
                               featureWeights = NULL, alpha = 0.5,
                               nfolds = 10, seed = 42) {
  if (!is.matrix(dataMatrix) || !is.numeric(dataMatrix)) {
    stop("dataMatrix must be a numeric matrix.")
  }
  if (is.null(rownames(dataMatrix))) {
    stop("dataMatrix must have row names corresponding to feature IDs.")
  }
  if (!is.numeric(responseVector)) {
    stop("responseVector must be numeric.")
  }
  if (length(responseVector) != ncol(dataMatrix)) {
    stop("Length of responseVector must match ncol(dataMatrix).")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha < 0 || alpha > 1) {
    stop("alpha must be a single numeric value between 0 and 1.")
  }
  if (!is.numeric(nfolds) || length(nfolds) != 1L || nfolds < 2 ||
      nfolds > ncol(dataMatrix)) {
    stop("nfolds must be between 2 and the number of samples.")
  }

  message("Step 1: Validating inputs and matching features...")
  validFeatures <- intersect(rownames(dataMatrix), candidateFeatures)
  if (length(validFeatures) == 0) {
    stop("No overlapping features between dataMatrix and candidateFeatures.")
  }
  message(sprintf("Found %d valid candidate features.", length(validFeatures)))

  if (is.null(featureWeights)) {
    message("No featureWeights provided. Training an unweighted model.")
    penalties <- rep(1, length(validFeatures))
  } else {
    if (!is.numeric(featureWeights) || is.null(names(featureWeights))) {
      stop("featureWeights must be a named numeric vector.")
    }
    missingWeights <- setdiff(validFeatures, names(featureWeights))
    if (length(missingWeights) > 0) {
      stop(sprintf("Missing weights for %d features.", length(missingWeights)))
    }

    penalties <- unname(featureWeights[validFeatures])
    if (any(!is.finite(penalties)) || any(penalties < 0)) {
      stop("featureWeights must contain finite non-negative values.")
    }
  }

  message("Step 2: Preparing data matrices for model training...")
  xTrain <- t(dataMatrix[validFeatures, , drop = FALSE])
  yTrain <- responseVector

  message("Step 3: Training Weighted Elastic Net model via cv.glmnet...")
  set.seed(seed)
  cvFit <- glmnet::cv.glmnet(
    x = xTrain,
    y = yTrain,
    family = "gaussian",
    alpha = alpha,
    standardize = TRUE,
    nfolds = nfolds,
    keep = TRUE,
    penalty.factor = penalties
  )

  message("Step 4: Extracting optimal model coefficients...")
  bestLambda <- cvFit$lambda.min
  coefMatrix <- as.matrix(stats::coef(cvFit, s = "lambda.min"))
  keepIdx <- which(coefMatrix[, 1] != 0)

  clockDf <- data.frame(
    Feature = rownames(coefMatrix)[keepIdx],
    Coefficient = unname(coefMatrix[keepIdx, 1]),
    stringsAsFactors = FALSE
  )

  message("Training complete.")
  list(
    cvModel = cvFit,
    bestLambda = bestLambda,
    clockCoefficients = clockDf
  )
}