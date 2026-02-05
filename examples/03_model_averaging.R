#' =============================================================================
#' GLmom Package - Example 3: Model Averaging for High Quantiles (MAGEV)
#' =============================================================================
#' This script demonstrates model averaging GEV estimation using ma.gev()
#' for robust high quantile estimation.
#'
#' Reference: Shin et al. (2026). Model averaging with mixed criteria for
#'            estimating high quantiles of extreme values: Application to
#'            heavy rainfall. SERA, 40(2), 47. doi:10.1007/s00477-025-03167-x
#' =============================================================================

# Load package
library(GLmom)

# Load example data
data(streamflow)
x <- streamflow$r1

cat("=== GLmom Example 3: Model Averaging (MAGEV) ===\n\n")
cat("Data: streamflow$r1 (n =", length(x), "annual maxima)\n\n")

# -----------------------------------------------------------------------------
# 1. Basic model averaging with default settings
# -----------------------------------------------------------------------------
cat("--- 1. Basic model averaging (likelihood weights) ---\n")
set.seed(123)
result <- ma.gev(x, quant = c(0.95, 0.99, 0.995), weight = "like1", B = 100)

cat("\nQuantile estimates comparison:\n")
comparison <- data.frame(
  Quantile = c("95%", "99%", "99.5%"),
  Return_Period = c("20-year", "100-year", "200-year"),
  MLE = round(result$qua.mle, 2),
  LME = round(result$qua.lme, 2),
  Model_Avg = round(result$qua.ma, 2)
)
print(comparison, row.names = FALSE)

cat("\nStandard errors (bootstrap, B=100):\n")
cat("  Fixed weights SE:", round(result$fixw.se.ma, 2), "\n")
cat("  Random weights SE:", round(result$ranw.se.ma, 2), "\n")

# -----------------------------------------------------------------------------
# 2. Compare different weighting schemes
# -----------------------------------------------------------------------------
cat("\n--- 2. Comparison of weighting schemes ---\n")
cat("Computing estimates for different weight methods...\n\n")

# Weight schemes to compare
weights <- c("like1", "gLd1", "med1", "cvt")
weight_names <- c("Likelihood (AIC)", "L-moment distance", "Median-based", "Conventional AIC")

results_99 <- data.frame(
  Weight = weight_names,
  Method = weights,
  Q99 = numeric(length(weights))
)

for (i in seq_along(weights)) {
  set.seed(123)
  res <- ma.gev(x, quant = 0.99, weight = weights[i], B = 50)
  results_99$Q99[i] <- res$qua.ma
}

cat("100-year return level by weighting scheme:\n")
print(data.frame(
  Weight_Method = results_99$Weight,
  Q99 = round(results_99$Q99, 2)
), row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. Trimming levels
# -----------------------------------------------------------------------------
cat("\n--- 3. Effect of trimming level ---\n")
cat("Trimming removes smallest observations from weight calculation.\n\n")

# Compare trimming levels for gLd weights
trims <- c("gLd0", "gLd1", "gLd2")
trim_levels <- c(0, 1, 2)

cat("99% quantile by trimming level (gLd weights):\n")
for (i in seq_along(trims)) {
  set.seed(123)
  res <- ma.gev(x, quant = 0.99, weight = trims[i], B = 50)
  cat("  trim =", trim_levels[i], ":", round(res$qua.ma, 2), "\n")
}

# -----------------------------------------------------------------------------
# 4. Model weights distribution
# -----------------------------------------------------------------------------
cat("\n--- 4. Model weights distribution ---\n")
set.seed(123)
result <- ma.gev(x, quant = 0.99, weight = "like1", B = 50)

cat("\nNumber of submodels (K):", length(result$w.ma), "\n")
cat("Weight distribution:\n")
cat("  Min:", round(min(result$w.ma), 4), "\n")
cat("  Max:", round(max(result$w.ma), 4), "\n")
cat("  Sum:", round(sum(result$w.ma), 4), "\n")

# Show top 5 weights with their xi values
if (length(result$pick_xi) == length(result$w.ma)) {
  weight_df <- data.frame(
    xi = result$pick_xi,
    weight = result$w.ma
  )
  weight_df <- weight_df[order(-weight_df$weight), ]
  cat("\nTop 5 submodels by weight:\n")
  print(head(data.frame(
    xi = round(weight_df$xi, 4),
    weight = round(weight_df$weight, 4)
  ), 5), row.names = FALSE)
}

# -----------------------------------------------------------------------------
# 5. Bayesian Model Averaging (BMA)
# -----------------------------------------------------------------------------
cat("\n--- 5. Bayesian Model Averaging ---\n")
set.seed(123)
result_bma <- ma.gev(x, quant = c(0.99, 0.995), weight = "like1",
                     bma = TRUE, pen = "norm", B = 50)

cat("\nBMA vs standard MA comparison (99% quantile):\n")
cat("  Standard MA:", round(result$qua.ma, 2), "\n")
cat("  BMA:", round(result_bma$qua.bma[1], 2), "\n")

# -----------------------------------------------------------------------------
# 6. Detailed output exploration
# -----------------------------------------------------------------------------
cat("\n--- 6. Available output components ---\n")
cat("\nMain outputs from ma.gev():\n")
cat("  $qua.ma      - Model-averaged quantile estimates\n")
cat("  $qua.mle    - MLE quantile estimates\n")
cat("  $qua.lme    - L-moment quantile estimates\n")
cat("  $w.ma       - Model weights\n")
cat("  $pick_xi    - Candidate shape parameters\n")
cat("  $fixw.se.ma - Bootstrap SE (fixed weights)\n")
cat("  $ranw.se.ma - Bootstrap SE (random weights)\n")
cat("  $qua.bma    - BMA quantiles (if bma=TRUE)\n")

# -----------------------------------------------------------------------------
# 7. Full example with all outputs
# -----------------------------------------------------------------------------
cat("\n--- 7. Complete analysis example ---\n")
set.seed(42)
full_result <- ma.gev(x,
                      quant = c(0.9, 0.95, 0.98, 0.99, 0.995),
                      weight = "gLd1",
                      B = 200,
                      bma = FALSE)

cat("\nFull return level table:\n")
full_table <- data.frame(
  Probability = c(0.9, 0.95, 0.98, 0.99, 0.995),
  Return_Period = c(10, 20, 50, 100, 200),
  MLE = round(full_result$qua.mle, 2),
  LME = round(full_result$qua.lme, 2),
  MA = round(full_result$qua.ma, 2),
  SE = round(full_result$fixw.se.ma, 2)
)
print(full_table, row.names = FALSE)

# -----------------------------------------------------------------------------
# 8. Coles-Dixon Penalized MLE (CD option)
# -----------------------------------------------------------------------------
cat("\n--- 8. Coles-Dixon Penalized MLE ---\n")
cat("CD=TRUE adds Coles-Dixon penalized MLE for shape parameter regularization.\n\n")

set.seed(123)
result_cd <- ma.gev(x, quant = c(0.99, 0.995), weight = "like1",
                    B = 50, CD = TRUE)

cat("Standard MLE:\n")
cat("  par:", round(result_cd$mle.hosking, 4), "\n")
cat("  99% quantile:", round(result_cd$qua.mle[1], 2), "\n")

cat("\nColes-Dixon penalized MLE:\n")
cat("  par:", round(result_cd$mle.CD, 4), "\n")
cat("  99% quantile:", round(result_cd$qua.CD[1], 2), "\n")

# -----------------------------------------------------------------------------
# 9. Restricted MLE (REMLE option)
# -----------------------------------------------------------------------------
cat("\n--- 9. Restricted MLE (REMLE) ---\n")
cat("remle=TRUE computes REMLE with mean constraint (stage 1) and mean+L-scale constraints (stage 2).\n\n")

set.seed(123)
result_remle <- ma.gev(x, quant = c(0.99, 0.995), weight = "like1",
                       B = 50, remle = TRUE)

cat("Standard MLE 99% quantile:", round(result_remle$qua.mle[1], 2), "\n")
cat("REMLE (mean constraint) 99% quantile:", round(result_remle$qua.remle1[1], 2), "\n")
cat("REMLE (mean+L-scale constraints) 99% quantile:", round(result_remle$qua.remle2[1], 2), "\n")

cat("\nREMLE parameters:\n")
cat("  remle1 (mean):", round(result_remle$remle1, 4), "\n")
cat("  remle2 (mean+L-scale):", round(result_remle$remle2, 4), "\n")

# -----------------------------------------------------------------------------
# 10. Combined CD and REMLE
# -----------------------------------------------------------------------------
cat("\n--- 10. Combined CD and REMLE estimation ---\n")
set.seed(123)
result_both <- ma.gev(x, quant = c(0.98, 0.99, 0.995), weight = "like1",
                      B = 50, CD = TRUE, remle = TRUE)

cat("\nQuantile comparison (all methods):\n")
comparison_all <- data.frame(
  Quantile = c("98%", "99%", "99.5%"),
  MLE = round(result_both$qua.mle, 2),
  CD = round(result_both$qua.CD, 2),
  REMLE_mean = round(result_both$qua.remle1, 2),
  REMLE_L1L2 = round(result_both$qua.remle2, 2),
  MA = round(result_both$qua.ma, 2)
)
print(comparison_all, row.names = FALSE)

# -----------------------------------------------------------------------------
# 11. Output quant field
# -----------------------------------------------------------------------------
cat("\n--- 11. Accessing requested quantiles ---\n")
cat("The 'quant' field returns the requested quantiles for convenience:\n")
cat("  result$quant:", result_both$quant, "\n")

cat("\n=== Example 3 completed ===\n")
