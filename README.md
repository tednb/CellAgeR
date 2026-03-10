# CellAgeR

CellAgeR provides enhancer-enriched, cell-type-specific epigenetic clocks for estimating DNA methylation age from beta value matrices. The package includes pre-trained clock coefficients, packaged example datasets, and utilities for computing both DNAm age and age acceleration.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("tednb/CellAgeR")
```

For local development from a source checkout:

```r
install.packages("devtools")
devtools::install(".", build_vignettes = TRUE)
```

## Quick Start

```r
library(CellAgeR)

clock_coefficients <- getClockCoefficients()
example_dir <- system.file("extdata", package = "CellAgeR")
load(file.path(example_dir, "Walker.Rdata"))

res <- predictCTSAge(
  betaMatrix = as.matrix(betaMat),
  clock = clock_coefficients$NeuronClock,
  age = phenoDf$Age
)

head(res$predictions)
head(res$ageAcceleration)
```

## License

CellAgeR is released under the MIT license. See `LICENSE` for details.