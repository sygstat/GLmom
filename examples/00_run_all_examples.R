#' =============================================================================
#' GLmom Package - Run All Examples
#' =============================================================================
#' This script runs all GLmom package examples sequentially.
#' Each example can also be run independently.
#'
#' Examples:
#'   01_stationary_gev.R     - Stationary GEV estimation with glme.gev()
#'   02_nonstationary_gev11.R - Non-stationary GEV11 with glme.gev11()
#'   03_model_averaging.R    - Model averaging with ma.gev()
#'   04_compatibility_functions.R - nsgev() and gado.prop_11()
#'   05_magev_diagnostics.R  - MAGEV diagnostic plots and new datasets
#'
#' Usage:
#'   source("00_run_all_examples.R")
#'
#' Or run individual examples:
#'   source("01_stationary_gev.R")
#' =============================================================================

cat("
================================================================================
                    GLmom Package - Complete Examples
================================================================================

This script runs all examples for the GLmom R package.
Package: GLmom - Generalized L-Moment Estimation for Extreme Value Distributions

References:
- Shin et al. (2025a). arXiv:2512.20385 (GLME)
- Shin et al. (2025b). JKSS, 54, 947-970 (NSGEV)
- Shin et al. (2026). SERA, 40(2), 47 (MAGEV)

================================================================================
")

# Check if GLmom is installed
if (!requireNamespace("GLmom", quietly = TRUE)) {
  stop("GLmom package is not installed. Please install it first:\n",
       "  remotes::install_github('sygstat/GLmom')")
}

# Get the directory of this script
script_dir <- dirname(sys.frame(1)$ofile)
if (is.null(script_dir) || script_dir == "") {
  script_dir <- getwd()
}

# List of example scripts
examples <- c(
  "01_stationary_gev.R",
  "02_nonstationary_gev11.R",
  "03_model_averaging.R",
  "04_compatibility_functions.R",
  "05_magev_diagnostics.R"
)

# Run each example
for (ex in examples) {
  cat("\n")
  cat(rep("=", 80), "\n", sep = "")
  cat("Running:", ex, "\n")
  cat(rep("=", 80), "\n", sep = "")

  ex_path <- file.path(script_dir, ex)
  if (file.exists(ex_path)) {
    tryCatch({
      source(ex_path, local = TRUE)
    }, error = function(e) {
      cat("\nError in", ex, ":", conditionMessage(e), "\n")
    })
  } else {
    cat("File not found:", ex_path, "\n")
  }

  cat("\n")
}

cat("
================================================================================
                         All Examples Completed
================================================================================

For more information:
- Package documentation: help(package = 'GLmom')
- GitHub: https://github.com/sygstat/GLmom
- Contact: syg.stat@etri.re.kr

================================================================================
")
