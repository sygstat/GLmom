#' =============================================================================
#' GLmom Package - Example 2: Non-stationary GEV11 Estimation
#' =============================================================================
#' This script demonstrates non-stationary GEV estimation using glme.gev11()
#' with time-varying location and scale parameters.
#'
#' Model: GEV11
#'   mu(t) = mu0 + mu1 * t
#'   sigma(t) = exp(sigma0 + sigma1 * t)
#'   xi = constant
#'
#' Reference: Shin et al. (2025b). Building nonstationary extreme value model
#'            using L-moments. Journal of the Korean Statistical Society, 54.
#' =============================================================================

# Load package
library(GLmom)

# Load example data: 53 years of river flow from Trehafod, Wales
data(Trehafod)
x <- Trehafod$r1
years <- Trehafod$Year

cat("=== GLmom Example 2: Non-stationary GEV11 Estimation ===\n\n")
cat("Data: Trehafod river flow (n =", length(x), "years)\n")
cat("Period:", min(years), "-", max(years), "\n")
cat("Summary:", "min =", round(min(x), 2), ", max =", round(max(x), 2),
    ", mean =", round(mean(x), 2), "\n\n")

# -----------------------------------------------------------------------------
# 1. Default estimation with adaptive beta penalty
# -----------------------------------------------------------------------------
cat("--- 1. GEV11 estimation with beta penalty (default) ---\n")
set.seed(123)  # For reproducibility
result <- glme.gev11(x, ntry = 10)

cat("\nGLME estimates (with penalty):\n")
cat("  mu0 =", round(result$para.glme[1], 4), "(baseline location)\n")
cat("  mu1 =", round(result$para.glme[2], 4), "(location trend per year)\n")
cat("  sigma0 =", round(result$para.glme[3], 4), "(log-scale baseline)\n")
cat("  sigma1 =", round(result$para.glme[4], 4), "(log-scale trend)\n")
cat("  xi =", round(result$para.glme[5], 4), "(shape parameter)\n")

cat("\nPure L-moment estimates (para.lme, Shin et al. 2025b method):\n")
cat("  mu0 =", round(result$para.lme[1], 4), "\n")
cat("  mu1 =", round(result$para.lme[2], 4), "\n")
cat("  sigma0 =", round(result$para.lme[3], 4), "\n")
cat("  sigma1 =", round(result$para.lme[4], 4), "\n")
cat("  xi =", round(result$para.lme[5], 4), "\n")

cat("\nGN16 method estimates (para.gado):\n")
cat("  mu0 =", round(result$para.gado[1], 4), "\n")
cat("  mu1 =", round(result$para.gado[2], 4), "\n")
cat("  sigma0 =", round(result$para.gado[3], 4), "\n")
cat("  sigma1 =", round(result$para.gado[4], 4), "\n")
cat("  xi =", round(result$para.gado[5], 4), "\n")

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
mu0 <- result$para.glme[1]
mu1 <- result$para.glme[2]
sigma0 <- result$para.glme[3]
sigma1 <- result$para.glme[4]
xi <- result$para.glme[5]

# 100-year return level at different time points
p <- 0.99  # 100-year return period
times <- c(1, 27, 53)  # start, middle, end
time_labels <- c("1968", "1994", "2020")

cat("\n100-year return levels over time:\n")
for (i in seq_along(times)) {
  t <- times[i]
  mu_t <- mu0 + mu1 * t
  sigma_t <- exp(sigma0 + sigma1 * t)

  # Compute quantile using GEV formula
  if (abs(xi) > 1e-6) {
    yp <- -log(p)
    q <- mu_t + sigma_t * (1 - yp^xi) / xi
  } else {
    q <- mu_t - sigma_t * log(-log(p))
  }
  cat("  Year", time_labels[i], "(t=", t, "):", round(q, 2), "\n")
}

# -----------------------------------------------------------------------------
# 6. Visualization (if running interactively)
# -----------------------------------------------------------------------------
cat("\n--- 6. Visualization ---\n")

# Plot data with fitted location trend
if (interactive()) {
  par(mfrow = c(1, 2))

  # Plot 1: Data with location trend
  plot(years, x, pch = 16, xlab = "Year", ylab = "Flow",
       main = "Trehafod River Flow with GEV11 Fit")
  t_seq <- 1:length(x)
  mu_fit <- mu0 + mu1 * t_seq
  lines(years, mu_fit, col = "blue", lwd = 2)
  legend("topleft", legend = c("Data", "mu(t)"),
         pch = c(16, NA), lty = c(NA, 1), col = c("black", "blue"))

  # Plot 2: Time-varying scale
  sigma_fit <- exp(sigma0 + sigma1 * t_seq)
  plot(years, sigma_fit, type = "l", col = "red", lwd = 2,
       xlab = "Year", ylab = "Scale parameter",
       main = "Time-varying Scale Parameter")

  par(mfrow = c(1, 1))
  cat("Plots displayed.\n")
} else {
  cat("Run interactively to see plots.\n")
}

cat("\n=== Example 2 completed ===\n")
