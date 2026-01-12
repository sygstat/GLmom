#' =============================================================================
#' GLmom Package - Example 5: MAGEV Diagnostic Visualization
#' =============================================================================
#' This script demonstrates the diagnostic plotting functions for MAGEV:
#' - magev.ksensplot(): K sensitivity analysis for optimal submodel selection
#' - magev.qqplot(): Q-Q diagnostic plot comparing estimation methods
#' - magev.rlplot(): Return level plot with confidence intervals
#'
#' Datasets used:
#' - bangkok: Annual maximum daily rainfall from Bangkok, Thailand
#' - haenam: Annual maximum daily rainfall from Haenam, South Korea
#'
#' Reference: Shin et al. (2025c). Model averaging with mixed criteria for
#'            estimating high quantiles of extreme values.
#'            arXiv:2505.21417
#' =============================================================================

# Load package
library(GLmom)

cat("=== GLmom Example 5: MAGEV Diagnostic Visualization ===\n\n")

# -----------------------------------------------------------------------------
# 1. Load Bangkok dataset
# -----------------------------------------------------------------------------
cat("--- 1. Bangkok Rainfall Data ---\n")
data(bangkok)
data_bkk <- bangkok[, 1]

cat("Data: Bangkok annual maximum daily rainfall\n")
cat("  Sample size: n =", length(data_bkk), "\n")
cat("  Range:", round(min(data_bkk), 1), "-", round(max(data_bkk), 1), "mm\n")
cat("  Mean:", round(mean(data_bkk), 1), "mm\n")
cat("  Source: Thai Meteorological Department (TMD)\n\n")

# -----------------------------------------------------------------------------
# 2. K Sensitivity Analysis (magev.ksensplot)
# -----------------------------------------------------------------------------
cat("--- 2. K Sensitivity Analysis ---\n")
cat("magev.ksensplot() helps determine the optimal number of submodels K.\n")
cat("It plots return level estimates, standard errors, and first-order\n")
cat("differences across K values from mink to maxk.\n\n")

if (interactive()) {
  # Run K sensitivity analysis
  cat("Running K sensitivity analysis (K = 4 to 20)...\n")
  optimal_k <- magev.ksensplot(data = data_bkk, q.cut = 0.6,
                                mink = 4, maxk = 20,
                                quant = c(0.99, 0.995))
  cat("\nOptimal K selected:", optimal_k, "\n")
  cat("(Optimal K minimizes oscillation in return level estimates)\n")
} else {
  cat("Run interactively to see K sensitivity plot.\n")
  cat("Example usage:\n")
  cat("  optimal_k <- magev.ksensplot(data_bkk, q.cut = 0.6,\n")
  cat("                                mink = 4, maxk = 20)\n")
  optimal_k <- 9  # Default for non-interactive
}

# -----------------------------------------------------------------------------
# 3. Q-Q Diagnostic Plot (magev.qqplot)
# -----------------------------------------------------------------------------
cat("\n--- 3. Q-Q Diagnostic Plot ---\n")
cat("magev.qqplot() creates a 2x2 panel of Q-Q plots comparing:\n")
cat("  - MLE (Maximum Likelihood Estimation)\n")
cat("  - LME (L-Moment Estimation)\n")
cat("  - Surrogate (Model-Averaged parameters)\n")
cat("  - REMLE (Restricted MLE with mean constraint)\n\n")

# Compute MAGEV estimates with remle option
cat("Computing MAGEV estimates with remle=TRUE...\n")
qq_probs <- c(seq(0.01, 0.99, by = 0.01), 0.995, 0.999)
set.seed(123)
zx_qq <- ma.gev(data = data_bkk, quant = qq_probs, weight = 'like1',
                numk = 9, varcom = FALSE, remle = TRUE, B = 50)

cat("Model-averaged 99% quantile:", round(zx_qq$zp.ma[which.min(abs(qq_probs - 0.99))], 2), "mm\n")
cat("REMLE (mean) 99% quantile:", round(zx_qq$qua.remle1[which.min(abs(qq_probs - 0.99))], 2), "mm\n\n")

if (interactive()) {
  cat("Generating Q-Q diagnostic plot...\n")
  magev.qqplot(data = data_bkk, zx = zx_qq)
  cat("Q-Q plot displayed.\n")
} else {
  cat("Run interactively to see Q-Q diagnostic plot.\n")
  cat("Example usage:\n")
  cat("  magev.qqplot(data = data_bkk, zx = zx_qq)\n")
}

# -----------------------------------------------------------------------------
# 4. Return Level Plot (magev.rlplot)
# -----------------------------------------------------------------------------
cat("\n--- 4. Return Level Plot ---\n")
cat("magev.rlplot() displays fitted return levels with 95% confidence intervals.\n")
cat("X-axis: Return period (log scale)\n")
cat("Y-axis: Return level (quantile)\n\n")

# Compute MAGEV estimates with variance components
cat("Computing MAGEV estimates with varcom=TRUE...\n")
ff <- c(seq(0.01, 0.09, by = 0.01), seq(0.1, 0.9, by = 0.1),
        0.93, 0.95, 0.98, 0.99, 0.995, 0.999)
set.seed(123)
zx_rl <- ma.gev(data = data_bkk, quant = ff, weight = 'like1',
                numk = 9, varcom = TRUE, B = 100)

cat("Surrogate GEV parameters:\n")
cat("  mu =", round(zx_rl$surr$par[1], 2), "\n")
cat("  sigma =", round(zx_rl$surr$par[2], 2), "\n")
cat("  xi =", round(zx_rl$surr$par[3], 4), "\n\n")

if (interactive()) {
  cat("Generating return level plot...\n")
  magev.rlplot(par = zx_rl$surr$par, se.vec = zx_rl$adj.se.ma, data = data_bkk)
  cat("Return level plot displayed.\n")
} else {
  cat("Run interactively to see return level plot.\n")
  cat("Example usage:\n")
  cat("  magev.rlplot(par = zx_rl$surr$par, se.vec = zx_rl$adj.se.ma, data = data_bkk)\n")
}

# -----------------------------------------------------------------------------
# 5. Haenam Dataset Example
# -----------------------------------------------------------------------------
cat("\n--- 5. Haenam Rainfall Data ---\n")
data(haenam)
data_haenam <- haenam[, 1]

cat("Data: Haenam annual maximum daily rainfall\n")
cat("  Sample size: n =", length(data_haenam), "\n")
cat("  Range:", round(min(data_haenam), 1), "-", round(max(data_haenam), 1), "mm\n")
cat("  Mean:", round(mean(data_haenam), 1), "mm\n")
cat("  Source: Korea Meteorological Administration (KMA)\n\n")

# Compute MAGEV for Haenam
cat("Computing MAGEV estimates for Haenam data...\n")
set.seed(456)
result_haenam <- ma.gev(data = data_haenam, quant = c(0.98, 0.99, 0.995),
                        weight = 'like1', numk = 10, B = 100,
                        CD = TRUE, remle = TRUE)

cat("\nReturn level estimates for Haenam:\n")
haenam_table <- data.frame(
  Return_Period = c("50-year", "100-year", "200-year"),
  Probability = c(0.98, 0.99, 0.995),
  MLE = round(result_haenam$qua.mle, 1),
  CD_MLE = round(result_haenam$qua.CD, 1),
  REMLE_mean = round(result_haenam$qua.remle1, 1),
  MA = round(result_haenam$zp.ma, 1)
)
print(haenam_table, row.names = FALSE)

# -----------------------------------------------------------------------------
# 6. Comparison: Bangkok vs Haenam
# -----------------------------------------------------------------------------
cat("\n--- 6. Comparison: Bangkok vs Haenam ---\n")

set.seed(789)
result_bkk <- ma.gev(data = data_bkk, quant = c(0.98, 0.99, 0.995),
                     weight = 'like1', numk = 9, B = 100)

cat("\n100-year return level comparison:\n")
cat("  Bangkok:", round(result_bkk$zp.ma[2], 1), "mm\n")
cat("  Haenam:", round(result_haenam$zp.ma[2], 1), "mm\n")

cat("\nShape parameter (xi) comparison:\n")
cat("  Bangkok MLE xi:", round(result_bkk$mle.hosking[3], 4), "\n")
cat("  Haenam MLE xi:", round(result_haenam$mle.hosking[3], 4), "\n")

# -----------------------------------------------------------------------------
# 7. Saving Plots to Files (Example)
# -----------------------------------------------------------------------------
cat("\n--- 7. Saving Plots to Files ---\n")
cat("To save diagnostic plots to files, use:\n\n")

cat("# K sensitivity plot\n")
cat("pdf('magev_ksensplot.pdf', width=8, height=6)\n")
cat("magev.ksensplot(data_bkk, quant=c(0.99, 0.995))\n")
cat("dev.off()\n\n")

cat("# Q-Q diagnostic plot\n")
cat("pdf('magev_qqplot.pdf', width=8, height=8)\n")
cat("magev.qqplot(data=data_bkk, zx=zx_qq)\n")
cat("dev.off()\n\n")

cat("# Return level plot\n")
cat("pdf('magev_rlplot.pdf', width=8, height=6)\n")
cat("magev.rlplot(par=zx_rl$surr$par, se.vec=zx_rl$adj.se.ma, data=data_bkk)\n")
cat("dev.off()\n")

# -----------------------------------------------------------------------------
# 8. Summary of Plotting Functions
# -----------------------------------------------------------------------------
cat("\n--- 8. Summary of Plotting Functions ---\n")

cat("
Function           Purpose                              Key Parameters
------------------ ------------------------------------ ----------------------
magev.ksensplot()  Select optimal K (submodel count)    data, mink, maxk, quant
magev.qqplot()     Compare estimation methods           data, zx (from ma.gev)
magev.rlplot()     Show return levels with 95% CI       par, se.vec, data

Usage tips:
- Use magev.ksensplot() first to determine optimal K
- Use zx from ma.gev(remle=TRUE) for magev.qqplot()
- Use zx from ma.gev(varcom=TRUE) for magev.rlplot()
")

cat("\n=== Example 5 completed ===\n")
