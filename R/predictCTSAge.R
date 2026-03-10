#' @title Predict Epigenetic Age and Age Acceleration
#' @description This function predicts epigenetic age based on a given DNA
#'     methylation beta matrix and a clock model. Optionally, if chronological
#'     age is provided, it calculates age acceleration residuals.
#'
#' @param betaMatrix A numeric matrix of beta values, with 
#'     CpGs in rows and samples in columns.
#' @param clock A data frame representing the epigenetic clock model. It must
#'     contain 'Feature' and 'Coefficient' columns. The intercept should be
#'     indicated by '(Intercept)' in the 'Feature' column.
#' @param age A numeric vector of chronological ages for the samples. The
#'     order must correspond to the sample order in `betaMatrix`. Defaults
#'     to NULL.
#'
#' @return A list containing the following components:
#' \item{predictions}{A data frame with the predicted DNAmAge for each sample.}
#' \item{ageAcceleration}{A data frame with Age Acceleration Residuals for
#'     each sample. NULL if `age` is not provided.}
#'
#' @importFrom stats lm residuals
#' @export
predictCTSAge <- function(betaMatrix, clock, age = NULL) {
  
  # --- 1. Process clock model ---
  message("Step 1: Processing clock model...")
  isIntercept <- clock$Feature == "(Intercept)"
  modelIntercept <- clock$Coefficient[isIntercept]
  
  if (length(modelIntercept) == 0) {
    warning("No intercept found in clock model. Using 0.")
    modelIntercept <- 0
  } else {
    modelIntercept <- as.numeric(modelIntercept[1])
  }
  
  # Extract all features by removing the intercept row
  featuresDf <- clock[!isIntercept, ]
  if (nrow(featuresDf) == 0) {
    stop("No valid features found in the clock model aside from intercept.")
  }
  
  # --- 2. Calculate DNAmAge ---
  message("Step 2: Calculating Epigenetic Age (DNAmAge)...")
  commonFeatures <- intersect(rownames(betaMatrix), featuresDf$Feature)
  
  if (length(commonFeatures) == 0) {
    stop("No overlapping features between betaMatrix and clock model.")
  }
  
  message(sprintf(
    "Matched %d / %d features.", 
    length(commonFeatures), nrow(featuresDf)
  ))
  
  featureCoefsSubset <- featuresDf[match(commonFeatures, featuresDf$Feature), ]
  betaSubset <- betaMatrix[commonFeatures, , drop = FALSE]
  
  # Matrix multiplication for speed and precision
  featureContribution <- as.numeric(featureCoefsSubset$Coefficient %*% betaSubset)
  dnamAge <- modelIntercept + featureContribution
  
  predictionsDf <- data.frame(
    SampleID = colnames(betaMatrix),
    DNAmAge = dnamAge,
    row.names = colnames(betaMatrix)
  )
  
  # --- 3. Calculate Age Acceleration ---
  ageAccelDf <- NULL
  
  if (!is.null(age)) {
    message("Step 3: Calculating age acceleration...")
    
    if (length(age) != ncol(betaMatrix)) {
      stop("Length of 'age' must match number of columns in 'betaMatrix'.")
    }
    
    tempDf <- data.frame(
      SampleID = colnames(betaMatrix),
      ChronologicalAge = age,
      DNAmAge = dnamAge
    )
    
    validIndices <- !is.na(tempDf$ChronologicalAge)
    cleanDf <- tempDf[validIndices, ]
    nDropped <- sum(!validIndices)
    
    if (nDropped > 0) {
      message(sprintf("Excluded %d samples missing age data.", nDropped))
    }
    
    if (nrow(cleanDf) < 3) {
      warning("Too few valid samples (< 3), skipping regression.")
    } else {
      # Initialize vector with NA_real_ for type safety
      accelTotalFull <- rep(NA_real_, nrow(predictionsDf))
      
      lmTotal <- lm(DNAmAge ~ ChronologicalAge, data = cleanDf)
      accelTotalFull[validIndices] <- residuals(lmTotal)
      
      ageAccelDf <- data.frame(
        AgeAcceleration = accelTotalFull,
        row.names = rownames(predictionsDf)
      )
    }
  } else {
    message("Step 3: 'age' not provided, skipping acceleration calculation.")
  }
  
  message("--- Analysis Complete ---")
  
  # Implicit return as per Bioconductor standard
  list(
    predictions = predictionsDf,
    ageAcceleration = ageAccelDf
  )
}