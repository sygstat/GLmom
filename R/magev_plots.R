# MAGEV Plotting Functions
#
# Reference: Shin, Y., Shin, Y., & Park, J. S. (2026). Model averaging with
# mixed criteria for estimating high quantiles of extreme values: Application
# to heavy rainfall. Stochastic Environmental Research and Risk Assessment,
# 40(2), 47. https://doi.org/10.1007/s00477-025-03167-x


#' K Sensitivity Plot for MAGEV
#'
#' @description
#' Plots return level estimates, standard errors, and first-order differences
#' across different numbers of candidate submodels K. This helps identify
#' stable regions where estimates converge and select an optimal K value.
#'
#' @param data A numeric vector of data to be fitted (e.g., annual maxima).
#' @param q.cut Quantile cutoff for determining stability (default 0.6).
#' @param mink Minimum number of candidate submodels to test (default 4).
#' @param maxk Maximum number of candidate submodels to test (default 20).
#' @param quant The probabilities for high quantile estimation.
#'   Default is c(0.99, 0.995).
#'
#' @details
#' The function computes MAGEV estimates for K ranging from \code{mink} to
#' \code{maxk}. For each K, it calculates:
#' \itemize{
#'   \item Return level estimates (black points)
#'   \item Normalized standard errors (blue line)
#'   \item First-order differences (red line with triangles)
#' }
#'
#' The optimal K is selected as the smallest value where both the normalized
#' standard error and first-order difference are below their respective
#' \code{q.cut} quantile cutoffs. The selected K is indicated by a purple
#' vertical line.
#'
#' @return The optimal K value (integer) selected by the algorithm.
#'
#' @references
#' Shin, Y., Shin, Y., & Park, J. S. (2026). Model averaging with mixed criteria
#' for estimating high quantiles of extreme values: Application to heavy rainfall.
#' \emph{Stochastic Environmental Research and Risk Assessment}, 40(2), 47.
#' \doi{10.1007/s00477-025-03167-x}
#'
#' @seealso \code{\link{ma.gev}} for the main model averaging function.
#'
#' @author Yonggwan Shin, Seokkap Ko, Jihong Park, Yire Shin, Jeong-Soo Park
#'
#' @examples
#' \donttest{
#' data(haenam)
#' optimal_k <- magev.ksensplot(haenam$X1)
#' print(optimal_k)
#' }
#'
#' @export
magev.ksensplot <- function(data = NULL, q.cut = 0.6, mink = 4,
                            maxk = 20, quant = c(0.99, 0.995)) {

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  totk <- maxk - mink + 1
  test <- list()
  numq <- length(quant)
  diff <- rt <- se <- matrix(NA, maxk, numq)

  qq <- quant
  id12 <- seq(1, numq)

  bma <- FALSE
  weight <- "like1"
  start <- "mle"

  for (k in mink:maxk) {

    test[[k]] <- ma.gev(data = data, quant = qq, weight = weight,
                        numk = k, B = 200, varcom = TRUE, trim = 0,
                        fig = FALSE, bma = bma, pen = "norm", remle = FALSE)

    if (bma != TRUE) {
      rt[k, id12] <- test[[k]]$qua.ma[id12]
      se[k, id12] <- test[[k]]$ranw.se.ma[id12]
    } else if (bma == TRUE) {
      rt[k, id12] <- test[[k]]$qua.bma[id12]
      se[k, id12] <- test[[k]]$pred.se.bma[id12]
    }

    message("k = ", k, ", zp.ma = ", paste(round(rt[k, id12], 2), collapse = ", "),
            ", se = ", paste(round(se[k, id12], 2), collapse = ", "))

  } # end for k

  tr.diff <- rep(NA, maxk)
  for (i in (mink + 1):(maxk - 1)) {
    diff[i, ] <- abs(rt[i + 1, ] - rt[i, ]) + abs(rt[i, ] - rt[i - 1, ])
  }

  nina <- which(!is.na(rt[, 1]))
  U <- (se[, 1] - min(se[, 1], na.rm = TRUE)) /
       (max(se[, 1], na.rm = TRUE) - min(se[, 1], na.rm = TRUE))
  tr.se <- U * (max(rt[, 1], na.rm = TRUE) - min(rt[, 1], na.rm = TRUE)) +
           min(rt[, 1], na.rm = TRUE)

  dina <- which(!is.na(diff[, 1]))

  Ub <- (diff[dina, 1] - min(diff[dina, 1], na.rm = TRUE)) /
        (max(diff[dina, 1], na.rm = TRUE) - min(diff[dina, 1], na.rm = TRUE))
  tr.diff[(mink + 1):(maxk - 1)] <- Ub * (max(rt[, 1], na.rm = TRUE) -
                                          min(rt[, 1], na.rm = TRUE)) +
                                    min(rt[, 1], na.rm = TRUE)

  diff.cut <- stats::quantile(tr.diff[dina], probs = q.cut, na.rm = TRUE)
  id <- which(tr.diff <= diff.cut)
  se.cut <- stats::quantile(tr.se[nina], probs = q.cut, na.rm = TRUE)

  tid <- which(tr.se <= se.cut)
  jd <- sort(intersect(id, tid))
  selTk <- min(jd)

  max.selTk <- min(maxk, selTk + 2)
  min.selTk <- max(mink, selTk - 2)
  iid <- intersect(jd, seq(min.selTk, max.selTk))
  star <- iid[which.min(se[iid, 1])]

  my <- min(rt[, 1], na.rm = TRUE) * 0.997
  My <- max(rt[, 1], na.rm = TRUE) * 1.003
  graphics::plot(seq(mink, maxk), rt[nina, 1], xlim = c(mink - 2, maxk),
                 ylim = c(my, My),
                 xlab = "K", ylab = "return level", pch = 16, cex = 1)

  graphics::axis(side = 1, at = seq(mink - 1, maxk, 1))

  graphics::par(new = TRUE)
  graphics::plot(seq(mink, maxk), tr.se[nina], xlim = c(mink - 2, maxk),
                 ylim = c(my, My),
                 type = "l", xlab = "K", ylab = "return level", col = "blue")
  graphics::points(seq(mink, maxk), tr.se[nina], col = "blue")

  graphics::par(new = TRUE)
  graphics::plot(seq(mink + 1, maxk - 1), tr.diff[(mink + 1):(maxk - 1)],
                 xlim = c(mink - 2, maxk), ylim = c(my, My),
                 type = "l", col = "red",
                 xlab = "K", ylab = "return level")
  graphics::points(seq(mink + 1, maxk - 1), tr.diff[(mink + 1):(maxk - 1)],
                   col = "red", pch = 2, cex = 0.8)

  graphics::abline(h = diff.cut, lty = 3, lwd = 2, col = "darkorange")
  graphics::abline(h = se.cut, lty = 3, lwd = 2, col = "darkgreen")
  graphics::abline(v = star, col = "purple")

  graphics::text(mink - 1.3, se.cut * 1.0012, "q.6(SE_K)", cex = 0.8, col = "darkgreen")
  graphics::text(mink - 1.5, diff.cut * 1.0012, "q.6(d_K)", col = "darkorange", cex = 0.8)

  return(star)
}


#' Q-Q Diagnostic Plot for MAGEV
#'
#' @description
#' Creates a 2x2 panel of Q-Q plots comparing observed vs. fitted quantiles
#' for different estimation methods: MLE, LME, surrogate (MA), and REMLE.
#'
#' @param data A numeric vector of observed data.
#' @param zx A list object returned by \code{\link{ma.gev}} with \code{remle = TRUE}.
#'
#' @details
#' The function creates four Q-Q plots:
#' \itemize{
#'   \item Upper left: MLE (Maximum Likelihood Estimation)
#'   \item Upper right: LME (L-moment Estimation)
#'   \item Lower left: Surrogate MA (Model Averaging surrogate)
#'   \item Lower right: REMLE (Restricted MLE, if available)
#' }
#'
#' Points close to the 45-degree diagonal line indicate good model fit.
#'
#' @return NULL. The function produces a plot as a side effect.
#'
#' @references
#' Shin, Y., Shin, Y., & Park, J. S. (2026). Model averaging with mixed criteria
#' for estimating high quantiles of extreme values: Application to heavy rainfall.
#' \emph{Stochastic Environmental Research and Risk Assessment}, 40(2), 47.
#' \doi{10.1007/s00477-025-03167-x}
#'
#' @seealso \code{\link{ma.gev}} for the main model averaging function,
#'   \code{\link{magev.rlplot}} for return level plots.
#'
#' @author Yonggwan Shin, Seokkap Ko, Jihong Park, Yire Shin, Jeong-Soo Park
#'
#' @examples
#' data(haenam)
#' qq <- c(seq(0.01, 0.99, by = 0.01), 0.995, 0.999)
#' zx <- ma.gev(haenam$X1, quant = qq, weight = 'like1',
#'              numk = 9, varcom = FALSE, remle = TRUE)
#' magev.qqplot(data = haenam$X1, zx = zx)
#'
#' @export
magev.qqplot <- function(data = NULL, zx = NULL) {

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  dat <- data
  qorg <- ((1:length(dat)) - 0.35) / (length(dat))

  graphics::par(mfrow = c(2, 2))
  minx <- min(dat) * 0.95

  maxp <- lmomco::quagev(max(qorg), lmomco::vec2par(zx$mle.hosking, 'gev'))
  maxb <- lmomco::quagev(max(qorg), lmomco::vec2par(zx$surr$par, 'gev'))
  maxx <- max(dat, maxp, maxb) * 1.2

  graphics::plot(lmomco::quagev(qorg, lmomco::vec2par(zx$mle.hosking, 'gev')),
                 sort(dat),
                 ylab = "Empirical", xlab = "model GEV", main = "mle",
                 xlim = c(minx, maxx), ylim = c(minx, maxx), pch = 16)
  graphics::abline(0, 1, col = 4)

  graphics::plot(lmomco::quagev(qorg, lmomco::vec2par(zx$lme, 'gev')),
                 sort(dat),
                 ylab = "Empirical", xlab = "model GEV", main = "lme",
                 xlim = c(minx, maxx), ylim = c(minx, maxx), pch = 16)
  graphics::abline(0, 1, col = 4)

  graphics::plot(lmomco::quagev(qorg, lmomco::vec2par(zx$surr$par, 'gev')),
                 sort(dat),
                 ylab = "Empirical", xlab = "model GEV", main = "surr.ma",
                 xlim = c(minx, maxx), ylim = c(minx, maxx), pch = 16)
  graphics::abline(0, 1, col = 4)

  if (!is.null(zx$remle1)) {
    graphics::plot(lmomco::quagev(qorg, lmomco::vec2par(zx$remle1, 'gev')),
                   sort(dat),
                   ylab = "Empirical", xlab = "model GEV", main = "remle1",
                   xlim = c(minx, maxx), ylim = c(minx, maxx), pch = 16)
    graphics::abline(0, 1, col = 4)
  }

  invisible(NULL)
}


#' Return Level Plot for MAGEV
#'
#' @description
#' Displays fitted return levels with 95% confidence intervals against
#' return period on a log scale.
#'
#' @param par A numeric vector of GEV parameters (mu, sigma, xi) in Hosking style.
#'   Typically from \code{zx$surr$par} where \code{zx} is the output of
#'   \code{\link{ma.gev}}.
#' @param se.vec A numeric vector of standard errors for the quantile estimates
#'   corresponding to the plotting positions. Typically from \code{zx$ranw.se.ma}.
#' @param data A numeric vector of observed data (annual maxima).
#'
#' @details
#' The plot shows:
#' \itemize{
#'   \item Black line: Fitted return level curve
#'   \item Blue lines: 95% confidence interval (mean +/- 1.96 * SE)
#'   \item Black points: Observed data at empirical return periods
#' }
#'
#' The x-axis (return period) is on a log scale, ranging from 0.1 to 900 years.
#'
#' @return NULL. The function produces a plot as a side effect.
#'
#' @references
#' Shin, Y., Shin, Y., & Park, J. S. (2026). Model averaging with mixed criteria
#' for estimating high quantiles of extreme values: Application to heavy rainfall.
#' \emph{Stochastic Environmental Research and Risk Assessment}, 40(2), 47.
#' \doi{10.1007/s00477-025-03167-x}
#'
#' @seealso \code{\link{ma.gev}} for the main model averaging function,
#'   \code{\link{magev.qqplot}} for Q-Q diagnostic plots.
#'
#' @author Yonggwan Shin, Seokkap Ko, Jihong Park, Yire Shin, Jeong-Soo Park
#'
#' @examples
#' data(haenam)
#' ff <- c(seq(0.01, 0.09, by = 0.01), 0.1, 0.2, 0.3, 0.4, 0.5,
#'         0.6, 0.7, 0.8, 0.9, 0.93, 0.95, 0.98, 0.99,
#'         0.993, 0.995, 0.998, 0.999)
#' zx <- ma.gev(haenam$X1, quant = ff, weight = 'like1',
#'              numk = 9, varcom = TRUE)
#' magev.rlplot(par = zx$surr$par, se.vec = zx$ranw.se.ma, data = haenam$X1)
#'
#' @export
magev.rlplot <- function(par = NULL, se.vec = NULL, data = NULL) {

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  a <- par  # a[3] is Hosking-style xi
  f <- c(seq(0.01, 0.09, by = 0.01), 0.1, 0.2, 0.3, 0.4, 0.5,
         0.6, 0.7, 0.8, 0.9, 0.93, 0.95, 0.98, 0.99,
         0.993, 0.995, 0.998, 0.999)
  q <- lmomco::quagev(f, lmomco::vec2par(a, 'gev'))

  v <- se.vec  # variance vector for f-quantiles
  maxq <- max(data) * 2.0
  ylima <- c(min(data, q), maxq)

  graphics::plot(-1 / log(f), q, log = "x", type = "n", xlim = c(0.1, 900),
                 ylim = ylima, xlab = "Return Period", ylab = "Return Level")

  qorg <- ((1:length(data)) - 0.35) / (length(data))

  graphics::lines(-1 / log(f), q)
  graphics::lines(-1 / log(f), q + 1.96 * v, col = 4)
  graphics::lines(-1 / log(f), q - 1.96 * v, col = 4)

  graphics::points(-1 / log(qorg), sort(data), pch = 16)

  invisible(NULL)
}
