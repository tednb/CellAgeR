#' @title Calculate Cell-Type Specific Score
#' @description Identify cell-type-specific features by combining empirical
#'   Bayes differential statistics with the minimum absolute difference in mean
#'   beta values between a target cell type and all remaining cell types.
#'
#' @param betaMatrix A numeric matrix of methylation beta values with features
#'   in rows and samples in columns.
#' @param cellTypes A character or factor vector describing the cell type for
#'   each sample in `betaMatrix`.
#' @param targetCellType A single string naming the cell type to evaluate.
#' @param convertToM Logical. If `TRUE`, beta values are converted to M-values
#'   before the limma model is fitted.
#' @param epsilon Small numeric value used to replace 0 and 1 beta values prior
#'   to M-value conversion.
#'
#' @return A data frame with one row per feature and the columns `Feature`,
#'   `logFC`, `P_Value`, `FDR`, and `CSS`, sorted by descending CSS and then
#'   increasing FDR.
#'
#' @examples
#' set.seed(1)
#' beta_mat <- matrix(
#'   runif(10 * 6, min = 0.05, max = 0.95),
#'   nrow = 10,
#'   dimnames = list(paste0("cg", 1:10), paste0("Sample", 1:6))
#' )
#' cell_types <- c("Neuron", "Neuron", "Neuron", "Astro", "Astro", "Astro")
#'
#' css_res <- calculateCSS(
#'   betaMatrix = beta_mat,
#'   cellTypes = cell_types,
#'   targetCellType = "Neuron"
#' )
#'
#' head(css_res)
#' @export
calculateCSS <- function(betaMatrix, cellTypes, targetCellType,
                         convertToM = TRUE, epsilon = 1e-6) {
  if (!is.matrix(betaMatrix) || !is.numeric(betaMatrix)) {
    stop("betaMatrix must be a numeric matrix.")
  }
  if (is.null(rownames(betaMatrix))) {
    stop("betaMatrix must have row names corresponding to feature IDs.")
  }

  cellTypes <- as.character(cellTypes)
  if (length(cellTypes) != ncol(betaMatrix)) {
    stop("Length of cellTypes must match the number of columns in betaMatrix.")
  }
  if (length(targetCellType) != 1L || !targetCellType %in% cellTypes) {
    stop("targetCellType must identify one cell type present in cellTypes.")
  }
  if (!is.logical(convertToM) || length(convertToM) != 1L || is.na(convertToM)) {
    stop("convertToM must be TRUE or FALSE.")
  }
  if (!is.numeric(epsilon) || length(epsilon) != 1L || is.na(epsilon) ||
      epsilon <= 0 || epsilon >= 0.5) {
    stop("epsilon must be a single numeric value between 0 and 0.5.")
  }

  uniqueCells <- unique(cellTypes)
  otherCells <- setdiff(uniqueCells, targetCellType)
  if (length(otherCells) == 0) {
    stop("Only one cell type found. At least two cell types are required.")
  }

  message("Calculating mean beta values per cell type...")
  meanBetaMat <- vapply(
    uniqueCells,
    FUN.VALUE = numeric(nrow(betaMatrix)),
    FUN = function(cell) {
      colsIdx <- which(cellTypes == cell)
      rowMeans(betaMatrix[, colsIdx, drop = FALSE], na.rm = TRUE)
    }
  )
  rownames(meanBetaMat) <- rownames(betaMatrix)

  message("Running empirical Bayes differential analysis (One-vs-Rest)...")
  fitMat <- betaMatrix
  if (convertToM) {
    fitMat[fitMat <= 0] <- epsilon
    fitMat[fitMat >= 1] <- 1 - epsilon
    fitMat <- log2(fitMat / (1 - fitMat))
  }

  cellFactor <- factor(cellTypes, levels = uniqueCells)
  design <- stats::model.matrix(~0 + cellFactor)
  colnames(design) <- levels(cellFactor)

  restExpr <- paste(otherCells, collapse = " + ")
  contrastFormula <- sprintf(
    "%s - (%s)/%d",
    targetCellType,
    restExpr,
    length(otherCells)
  )
  contMatrix <- limma::makeContrasts(contrasts = contrastFormula, levels = design)

  fit <- limma::lmFit(fitMat, design)
  fit2 <- limma::contrasts.fit(fit, contMatrix)
  fit2 <- limma::eBayes(fit2)
  limmaRes <- limma::topTable(fit2, coef = 1, number = Inf, sort.by = "none")

  message(sprintf("Calculating CSS for %s...", targetCellType))
  targetMeans <- meanBetaMat[, targetCellType]
  otherMeans <- meanBetaMat[, otherCells, drop = FALSE]
  cssScores <- matrixStats::rowMins(abs(targetMeans - otherMeans), na.rm = TRUE)

  finalDf <- data.frame(
    Feature = rownames(limmaRes),
    logFC = limmaRes$logFC,
    P_Value = limmaRes$P.Value,
    FDR = limmaRes$adj.P.Val,
    CSS = cssScores,
    stringsAsFactors = FALSE
  )

  finalDf <- finalDf[order(-finalDf$CSS, finalDf$FDR), ]
  rownames(finalDf) <- NULL

  message("CSS calculation complete.")
  finalDf
}