#' =============================================================================
#' GLmom Package - Example 4: Compatibility Functions
#' =============================================================================
#' This script demonstrates the compatibility functions nsgev() and
#' gado.prop_11() that follow the methodology in Shin et al. (2025b).
#'
#' These functions provide:
#' - nsgev(): Simple interface returning only the proposed L-moment estimates
#' - gado.prop_11(): Comprehensive output comparing multiple estimation methods
#'
#' Reference: Shin et al. (2025b). Building nonstationary extreme value model
#'            using L-moments. Journal of the Korean Statistical Society, 54.
#' =============================================================================

# Load package
library(GLmom)

# Load example data
data(Trehafod)
x <- Trehafod$r1

cat("=== GLmom Example 4: Compatibility Functions ===\n\n")
cat("Data: Trehafod river flow (n =", length(x), "years)\n\n")

# -----------------------------------------------------------------------------
# 1. nsgev() - Simple interface
# -----------------------------------------------------------------------------
cat("--- 1. nsgev() - Simple L-moment estimation ---\n")
cat("This function returns only the proposed L-moment estimates (no penalty).\n\n")

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
cat("\n--- 2. gado.prop_11() - Comprehensive comparison ---\n")
cat("This function returns estimates from three different methods.\n\n")

set.seed(123)
result2 <- gado.prop_11(x, ntry = 10)

cat("Method comparison:\n\n")

# Create comparison table
methods <- c("Proposed (L-moment)", "GN16", "WLS (strup.final)")
estimates <- rbind(
  result2$para.prop,
  result2$para.gado,
  result2$strup.final
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

cat("\nStationary WLS (strup.sta) - initial values:\n")
cat("  ", round(result2$strup.sta, 4), "\n")

cat("\nOriginal WLS (strup.org):\n")
cat("  ", round(result2$strup.org, 4), "\n")

cat("\nFinal WLS (strup.final):\n")
cat("  ", round(result2$strup.final, 4), "\n")

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
cat("\nglme.gev11(pen='no')$para.jkss:\n  ", round(result_glme$para.jkss, 4), "\n")

# Check if they match
if (all(abs(result1$para.prop - result_glme$para.jkss) < 1e-6)) {
  cat("\nThe results match exactly.\n")
} else {
  cat("\nNote: Results may differ slightly due to random initialization.\n")
}

# -----------------------------------------------------------------------------
# 5. When to use which function
# -----------------------------------------------------------------------------
cat("\n--- 5. Function selection guide ---\n")

cat("
Use nsgev() when:
  - You want the simplest interface

  - You only need the proposed L-moment estimates
  - You're following Shin et al. (2025b) methodology exactly

Use gado.prop_11() when:
  - You want to compare multiple estimation methods
  - You need WLS and GN16 estimates for diagnostics
  - You're conducting a methodological comparison study

Use glme.gev11() when:
  - You want penalized estimates (GLME)
  - You need full control over penalty functions
  - You want both penalized and unpenalized estimates
")

# -----------------------------------------------------------------------------
# 6. Example with PhliuAgromet data
# -----------------------------------------------------------------------------
cat("\n--- 6. Example with PhliuAgromet data ---\n")
data(PhliuAgromet)

# Check if precipitation data is available
if ("prec" %in% names(PhliuAgromet)) {
  y <- PhliuAgromet$prec
  y <- y[!is.na(y)]

  if (length(y) >= 20) {
    cat("Data: PhliuAgromet precipitation (n =", length(y), ")\n\n")

    set.seed(456)
    result_philu <- nsgev(y, ntry = 5)

    cat("Proposed estimates:\n")
    cat("  mu0 =", round(result_philu$para.prop[1], 4), "\n")
    cat("  mu1 =", round(result_philu$para.prop[2], 4), "\n")
    cat("  xi =", round(result_philu$para.prop[5], 4), "\n")
  } else {
    cat("Insufficient data in PhliuAgromet for this example.\n")
  }
} else {
  cat("PhliuAgromet$prec not available. Skipping this example.\n")
}

cat("\n=== Example 4 completed ===\n")
