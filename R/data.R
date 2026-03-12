#' Cell-Type Specific DNA Methylation Signatures
#'
#' A list containing named numeric vectors of Cell-Type Specific  
#' scores (CSS) for various cell types. These CpG features were selected based on 
#' a baseline threshold (CSS > 0.1) to provide a reference for 
#' downstream analyses.
#'
#' @format A list of 8 named numeric vectors of CSS scores, with Illumina CpG IDs
#' as names:
#' \describe{
#'   \item{cd4t}{CD4+ T cells}
#'   \item{mono}{Monocytes}
#'   \item{colonEpi}{Colon epithelium}
#'   \item{lungEpi}{Lung epithelium}
#'   \item{keratinocyte}{Keratinocytes}
#'   \item{hepatocyte}{Hepatocytes}
#'   \item{neuron}{Neurons}
#'   \item{oligo}{Oligodendrocytes}
#' }
#' 
#' @source 
#' The CSS scores were calculated using the `calculateCSS` function within this 
#' package, based on highly purified cell-type DNA methylation reference datasets 
#' provided by Loyfer et al. and Bell et al.
#' 
#' @usage data(cssSignatures)
"cssSignatures"