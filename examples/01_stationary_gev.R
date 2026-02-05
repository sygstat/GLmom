#' =============================================================================
#' GLmom Package - Example 1: Stationary GEV Estimation
#' =============================================================================
#' This script demonstrates stationary GEV parameter estimation using glme.gev()
#' with various penalty functions.
#'
#' Reference: Shin et al. (2025a). Generalized method of L-moment estimation
#'            for stationary and nonstationary extreme value models.
#'            arXiv:2512.20385
#' =============================================================================

# Load package
library(GLmom)

# Load example data: 50 annual maximum streamflow values
data(streamflow)
x <- streamflow$r1

cat("=== GLmom Example 1: Stationary GEV Estimation ===\n\n")
cat("Data: streamflow$r1 (n =", length(x), "annual maxima)\n")
cat("Summary:", "min =", round(min(x), 2), ", max =", round(max(x), 2),
    ", mean =", round(mean(x), 2), "\n\n")

# -----------------------------------------------------------------------------
# 1. Default estimation with adaptive beta penalty
# -----------------------------------------------------------------------------
cat("--- 1. Default estimation (adaptive beta penalty) ---\n")
result <- glme.gev(x)

cat("GLME estimates (with penalty):\n")
cat("  mu =", round(result$para.glme[1], 4), "\n")
cat("  sigma =", round(result$para.glme[2], 4), "\n")
cat("  xi =", round(result$para.glme[3], 4), "\n\n")

cat("L-moment estimates (no penalty):\n")
cat("  mu =", round(result$para.lme[1], 4), "\n")
cat("  sigma =", round(result$para.lme[2], 4), "\n")
cat("  xi =", round(result$para.lme[3], 4), "\n\n")

cat("Negative log-likelihood:", round(result$nllh.glme, 4), "\n\n")

# -----------------------------------------------------------------------------
# 2. Compare different penalty functions
# -----------------------------------------------------------------------------
cat("--- 2. Comparison of penalty functions ---\n")
penalties <- c("beta", "norm", "ms", "park", "cannon", "cd", "no")

results_df <- data.frame(
  penalty = penalties,
  mu = numeric(length(penalties)),
  sigma = numeric(length(penalties)),
  xi = numeric(length(penalties))
)

for (i in seq_along(penalties)) {
  res <- glme.gev(x, pen = penalties[i])
  results_df$mu[i] <- res$para.glme[1]
  results_df$sigma[i] <- res$para.glme[2]
  results_df$xi[i] <- res$para.glme[3]
}

cat("\nShape parameter (xi) estimates by penalty function:\n")
print(data.frame(
  penalty = results_df$penalty,
  xi = round(results_df$xi, 4)
), row.names = FALSE)

cat("\nFull parameter estimates:\n")
print(data.frame(
  penalty = results_df$penalty,
  mu = round(results_df$mu, 4),
  sigma = round(results_df$sigma, 4),
  xi = round(results_df$xi, 4)
), row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. Custom hyperparameters for beta penalty
# -----------------------------------------------------------------------------
cat("\n--- 3. Beta penalty with custom hyperparameters ---\n")

# Using pen.choice presets (1-6)
cat("Preset hyperparameter choices for beta penalty:\n")
for (choice in 1:6) {
  res <- glme.gev(x, pen = "beta", pen.choice = choice)
  cat("  pen.choice =", choice, ": xi =", round(res$para.glme[3], 4), "\n")
}

# Using custom p, c1, c2 values
cat("\nCustom hyperparameters (p=6, c1=20, c2=7):\n")
res_custom <- glme.gev(x, pen = "beta", p = 6, c1 = 20, c2 = 7)
cat("  xi =", round(res_custom$para.glme[3], 4), "\n")

# -----------------------------------------------------------------------------
# 4. Normal penalty with custom mean and std
# -----------------------------------------------------------------------------
cat("\n--- 4. Normal penalty with custom parameters ---\n")

# Preset choices for norm penalty (1-4)
cat("Preset hyperparameter choices for norm penalty:\n")
for (choice in 1:4) {
  res <- glme.gev(x, pen = "norm", pen.choice = choice)
  cat("  pen.choice =", choice, ": xi =", round(res$para.glme[3], 4), "\n")
}

# Custom normal prior
cat("\nCustom normal prior (mu=-0.5, std=0.2):\n")
res_norm <- glme.gev(x, pen = "norm", mu = -0.5, std = 0.2)
cat("  xi =", round(res_norm$para.glme[3], 4), "\n")

# -----------------------------------------------------------------------------
# 5. Quantile estimation
# -----------------------------------------------------------------------------
cat("\n--- 5. Quantile estimation ---\n")
result <- glme.gev(x, pen = "beta")

# Extract parameters
mu <- result$para.glme[1]
sigma <- result$para.glme[2]
xi <- result$para.glme[3]

# Compute return levels using lmomco
# lmomco GEV parameterization: xi=location, alpha=scale, kappa=shape
para <- lmomco::vec2par(c(mu, sigma, xi), type = "gev")
probs <- c(0.9, 0.95, 0.99, 0.995)

cat("Return levels (GLME estimates with beta penalty):\n")
for (p in probs) {
  q <- lmomco::quagev(p, para)
  T <- 1 / (1 - p)
  cat(sprintf("  %.0f-year return level (p=%.3f): %.2f\n", T, p, q))
}

cat("\n=== Example 1 completed ===\n")
