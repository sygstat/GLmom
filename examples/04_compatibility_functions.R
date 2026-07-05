#' =============================================================================
#' GLmom Package - Example 4: Compatibility Functions
#' =============================================================================
#' This script demonstrates the compatibility functions nsgev() and
#' gado.prop_11() that follow the methodology in Shin et al. (2025b).
#'
#' Since v2.0.0, gado.prop_11() is deprecated: the proposed method is
#' available as lme.gev11(), and the auxiliary methods as strup.gev11()
#' and GN16.gev11(). The wrappers below remain for v1.x compatibility.
#'
#' These functions provide:
#' - lme.gev11(): Recommended interface for the proposed L-moment method
#' - nsgev(): Simple wrapper returning only the proposed estimates
#' - gado.prop_11(): Comprehensive output comparing multiple methods (deprecated)
#'
#' Reference: Shin et al. (2025b). Building nonstationary extreme value model
#'            using L-moments. Journal of the Korean Statistical Society, 54.
#' =============================================================================

# Load package
library(GLmom)

# Load example data
data(PhliuAgromet)
x <- PhliuAgromet$prec

cat("=== GLmom Example 4: Compatibility Functions ===\n\n")
cat("Data: Phliu annual maximum precipitation (n =", length(x), "years)\n\n")

# -----------------------------------------------------------------------------
# 1. nsgev() - Simple interface
# -----------------------------------------------------------------------------
cat("--- 1. lme.gev11() / nsgev() - L-moment estimation ---\n")
cat("These return the proposed L-moment estimates (no penalty).\n\n")

# Recommended v2.0.0 interface:
set.seed(123)
result0 <- lme.gev11(x, ntry = 10)
cat("lme.gev11() estimates:\n  ", round(result0$lme.gev11, 4), "\n\n")

# v1.x-compatible wrapper (same method):
set.seed(123)
result1 <- nsgev(x, ntry = 10)

cat("Proposed L-moment estimates (para.prop):\n")
cat("  mu0 =", round(result1$para.prop[1], 4), "\n")
cat("  mu1 =", round(result1$para.prop[2], 4), "\n")
cat("  sigma0 =", round(result1$para.prop[3], 4), "\n")
cat("  sigma1 =", round(result1$para.prop[4], 4), "\n")
cat("  xi =", round(result1$para.prop[5], 4), "\n")

cat("\nPrecision measure:", round(result1$precis, 6), "\n")

# -----------------------------------------------------------------------------
# 2. gado.prop_11() - Comprehensive output
# -----------------------------------------------------------------------------
cat("\n--- 2. gado.prop_11() - Comprehensive comparison (deprecated) ---\n")
cat("This function returns estimates from three different methods.\n")
cat("(A deprecation warning is expected; use lme.gev11() going forward.)\n\n")

set.seed(123)
result2 <- suppressWarnings(gado.prop_11(x, ntry = 10))

cat("Method comparison:\n\n")

# Create comparison table
methods <- c("Proposed (L-moment)", "GN16", "WLS (para.wls)")
estimates <- rbind(
  result2$para.prop,
  result2$para.gado,
  result2$para.wls
)

comparison <- data.frame(
  Method = methods,
  mu0 = round(estimates[, 1], 4),
  mu1 = round(estimates[, 2], 4),
  sigma0 = round(estimates[, 3], 4),
  sigma1 = round(estimates[, 4], 4),
  xi = round(estimates[, 5], 4)
)
print(comparison, row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. Additional outputs from gado.prop_11()
# -----------------------------------------------------------------------------
cat("\n--- 3. Additional WLS estimates ---\n")

cat("\nOriginal WLS (strup.org):\n")
cat("  ", round(result2$strup.org, 4), "\n")

cat("\nFinal WLS (para.wls):\n")
cat("  ", round(result2$para.wls, 4), "\n")

cat("\nStationary L-moments (lme.sta) - (mu, sigma, xi):\n")
cat("  mu =", round(result2$lme.sta[1], 4), "\n")
cat("  sigma =", round(result2$lme.sta[2], 4), "\n")
cat("  xi =", round(result2$lme.sta[3], 4), "\n")

# -----------------------------------------------------------------------------
# 4. Relationship to glme.gev11()
# -----------------------------------------------------------------------------
cat("\n--- 4. Relationship to glme.gev11() ---\n")
cat("nsgev() and gado.prop_11() are wrappers around glme.gev11(pen='no').\n\n")

set.seed(123)
result_glme <- glme.gev11(x, ntry = 10, pen = "no")

cat("Comparison of outputs:\n")
cat("\nnsgev()$para.prop:\n  ", round(result1$para.prop, 4), "\n")
cat("\nglme.gev11(pen='no')$para.lme:\n  ", round(result_glme$para.lme, 4), "\n")

# Check if they match
if (all(abs(result1$para.prop - result_glme$para.lme) < 1e-6)) {
  cat("\nThe results match exactly.\n")
} else {
  cat("\nNote: Results may differ slightly due to random initialization.\n")
}

# -----------------------------------------------------------------------------
# 5. When to use which function
# -----------------------------------------------------------------------------
cat("\n--- 5. Function selection guide ---\n")

cat("
Use lme.gev11() when:
  - You want the proposed L-moment method of Shin et al. (2025b)
  - (nsgev() is an equivalent v1.x-compatible wrapper)

Use strup.gev11() / GN16.gev11() when:
  - You need the WLS or GN16 estimates for diagnostics
  - You're conducting a methodological comparison study
  - (gado.prop_11() bundles all three but is deprecated)

Use glme.gev11() when:
  - You want penalized estimates (GLME)
  - You need full control over penalty functions
")

# -----------------------------------------------------------------------------
# 6. Example with haenam data
# -----------------------------------------------------------------------------
cat("\n--- 6. Example with haenam data ---\n")
data(haenam)
y <- haenam$X1

cat("Data: Haenam annual maximum rainfall (n =", length(y), ")\n\n")

set.seed(456)
result_haenam <- nsgev(y, ntry = 5)

cat("Proposed estimates:\n")
cat("  mu0 =", round(result_haenam$para.prop[1], 4), "\n")
cat("  mu1 =", round(result_haenam$para.prop[2], 4), "\n")
cat("  xi =", round(result_haenam$para.prop[5], 4), "\n")

cat("\n=== Example 4 completed ===\n")
