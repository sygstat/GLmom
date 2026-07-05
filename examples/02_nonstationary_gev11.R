#' =============================================================================
#' GLmom Package - Example 2: Non-stationary GEV11 Estimation
#' =============================================================================
#' This script demonstrates non-stationary GEV estimation using glme.gev11()
#' with time-varying location and scale parameters, together with the
#' component methods lme.gev11(), strup.gev11(), and GN16.gev11()
#' (new in v2.0.0).
#'
#' Model: GEV11
#'   mu(t) = mu0 + mu1 * t
#'   sigma(t) = exp(sigma0 + sigma1 * t)
#'   xi = constant
#'
#' References:
#'   Shin et al. (2025a). Generalized method of L-moment estimation for
#'     stationary and nonstationary extreme value models. arXiv:2512.20385.
#'   Shin et al. (2025b). Building nonstationary extreme value model
#'     using L-moments. Journal of the Korean Statistical Society, 54.
#' =============================================================================

# Load package
library(GLmom)

# Load example data: 40 years of annual maximum daily precipitation at the
# Phliu Agrometeorological Station, Thailand (used in Shin et al. 2025a)
data(PhliuAgromet)
x <- PhliuAgromet$prec
years <- PhliuAgromet$year

cat("=== GLmom Example 2: Non-stationary GEV11 Estimation ===\n\n")
cat("Data: Phliu annual maximum precipitation (n =", length(x), "years)\n")
cat("Period:", min(years), "-", max(years), "\n")
cat("Summary:", "min =", round(min(x), 2), ", max =", round(max(x), 2),
    ", mean =", round(mean(x), 2), "\n\n")

# -----------------------------------------------------------------------------
# 1. Default estimation with adaptive beta penalty
# -----------------------------------------------------------------------------
cat("--- 1. GEV11 estimation with beta penalty (default) ---\n")
set.seed(123)  # For reproducibility
result <- glme.gev11(x, ntry = 10)

# v2.0.0: classed object ("glme11") -- print() gives a concise summary,
# plot() draws a Gumbel-scale Q-Q diagnostic (see section 7)
print(result)

cat("\nGLME estimates (with penalty):\n")
cat("  mu0 =", round(result$para.glme[1], 4), "(baseline location)\n")
cat("  mu1 =", round(result$para.glme[2], 4), "(location trend per year)\n")
cat("  sigma0 =", round(result$para.glme[3], 4), "(log-scale baseline)\n")
cat("  sigma1 =", round(result$para.glme[4], 4), "(log-scale trend)\n")
cat("  xi =", round(result$para.glme[5], 4), "(shape parameter)\n")

# Pure L-moment estimates: use lme.gev11() (v2.0.0; formerly gado.prop_11)
set.seed(123)
lme_res <- lme.gev11(x, ntry = 10)
cat("\nPure L-moment estimates (lme.gev11, Shin et al. 2025b method):\n")
cat("  mu0 =", round(lme_res$lme.gev11[1], 4), "\n")
cat("  mu1 =", round(lme_res$lme.gev11[2], 4), "\n")
cat("  sigma0 =", round(lme_res$lme.gev11[3], 4), "\n")
cat("  sigma1 =", round(lme_res$lme.gev11[4], 4), "\n")
cat("  xi =", round(lme_res$lme.gev11[5], 4), "\n")

cat("\nGN16 method estimates (para.gado):\n")
cat("  mu0 =", round(result$para.gado[1], 4), "\n")
cat("  mu1 =", round(result$para.gado[2], 4), "\n")
cat("  sigma0 =", round(result$para.gado[3], 4), "\n")
cat("  sigma1 =", round(result$para.gado[4], 4), "\n")
cat("  xi =", round(result$para.gado[5], 4), "\n")

cat("\nWLS estimates (para.wls, Strupczewski method):\n")
cat("  ", paste(round(result$para.wls, 4), collapse = "  "), "\n")

# -----------------------------------------------------------------------------
# 2. Interpretation of trend parameters
# -----------------------------------------------------------------------------
cat("\n--- 2. Interpretation of trend parameters ---\n")
mu0 <- result$para.glme[1]
mu1 <- result$para.glme[2]
sigma0 <- result$para.glme[3]
sigma1 <- result$para.glme[4]

# Time-varying parameters at start and end of period
t_start <- 1
t_end <- length(x)

mu_start <- mu0 + mu1 * t_start
mu_end <- mu0 + mu1 * t_end
sigma_start <- exp(sigma0 + sigma1 * t_start)
sigma_end <- exp(sigma0 + sigma1 * t_end)

cat("\nLocation parameter mu(t) = mu0 + mu1*t:\n")
cat("  At start (t=1):", round(mu_start, 2), "\n")
cat("  At end (t=", t_end, "):", round(mu_end, 2), "\n")
cat("  Total change:", round(mu_end - mu_start, 2), "\n")

cat("\nScale parameter sigma(t) = exp(sigma0 + sigma1*t):\n")
cat("  At start (t=1):", round(sigma_start, 2), "\n")
cat("  At end (t=", t_end, "):", round(sigma_end, 2), "\n")
cat("  Relative change:", round((sigma_end / sigma_start - 1) * 100, 1), "%\n")

# -----------------------------------------------------------------------------
# 3. Compare penalty functions for non-stationary model
# -----------------------------------------------------------------------------
cat("\n--- 3. Shape parameter under different penalties ---\n")
penalties <- c("beta", "norm", "ms", "park", "cannon", "cd", "no")

cat("\nShape parameter (xi) by penalty function:\n")
for (pen in penalties) {
  set.seed(123)
  res <- glme.gev11(x, ntry = 5, pen = pen)
  cat("  ", pen, ": xi =", round(res$para.glme[5], 4), "\n")
}

cat("\nShape parameter (xi) by beta penalty preset (pen.choice 1-6):\n")
for (ch in 1:6) {
  set.seed(123)
  res <- glme.gev11(x, ntry = 5, pen = "beta", pen.choice = ch)
  cat("  pen.choice =", ch, ": xi =", round(res$para.glme[5], 4), "\n")
}

# -----------------------------------------------------------------------------
# 4. Effect of ntry (multi-start optimization)
# -----------------------------------------------------------------------------
cat("\n--- 4. Effect of ntry on estimation stability ---\n")
cat("Running with different ntry values...\n")

ntry_values <- c(5, 10, 20, 30)
for (nt in ntry_values) {
  set.seed(42)
  res <- glme.gev11(x, ntry = nt, pen = "beta")
  cat("  ntry =", sprintf("%2d", nt), ": xi =", round(res$para.glme[5], 4), "\n")
}

# -----------------------------------------------------------------------------
# 5. Time-varying return levels
# -----------------------------------------------------------------------------
cat("\n--- 5. Time-varying return levels ---\n")
xi <- result$para.glme[5]

# 100-year return level at different time points (via quagev.NS)
q100 <- quagev.NS(f = 0.99, para = result$para.glme,
                  nsample = length(x), model = "gev11")

cat("\n100-year return levels over time:\n")
for (t in c(1, round(length(x)/2), length(x))) {
  cat("  Year", years[t], "(t=", t, "):", round(q100[t], 2), "\n")
}

# -----------------------------------------------------------------------------
# 6. Random sample generation from the fitted model (new in v2.0.0)
# -----------------------------------------------------------------------------
cat("\n--- 6. Random sample from the fitted GEV11 model ---\n")
set.seed(1)
sim <- ran.gev_all(length(x), para = result$para.glme, model = "gev11")
cat("Simulated series: mean =", round(mean(sim), 2),
    ", max =", round(max(sim), 2), "\n")

# -----------------------------------------------------------------------------
# 7. Visualization (if running interactively)
# -----------------------------------------------------------------------------
cat("\n--- 7. Visualization ---\n")

# Plot data with fitted location trend
if (interactive()) {
  par(mfrow = c(1, 2))

  # Plot 1: Data with location trend
  plot(years, x, pch = 16, xlab = "Year", ylab = "Precipitation (mm)",
       main = "Phliu Annual Maximum Precipitation with GEV11 Fit")
  t_seq <- 1:length(x)
  mu_fit <- mu0 + mu1 * t_seq
  lines(years, mu_fit, col = "blue", lwd = 2)
  legend("topleft", legend = c("Data", "mu(t)"),
         pch = c(16, NA), lty = c(NA, 1), col = c("black", "blue"))

  # Plot 2: Gumbel-scale Q-Q diagnostic via the plot method (v2.0.0)
  plot(result)

  par(mfrow = c(1, 1))
  cat("Plots displayed.\n")
} else {
  cat("Run interactively to see plots.\n")
}

cat("\n=== Example 2 completed ===\n")
