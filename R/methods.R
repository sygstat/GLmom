# S3 methods (print, summary, plot) for the fitted objects returned by
# glme.gev() ["glme"], glme.gev11() ["glme11"], lme.gev11() ["lme11"],
# and ma.gev() ["magev"]. All methods return their first argument invisibly.

# one-line description of the penalty settings of a fit -------------------
penalty.label <- function(x) {
  pen <- x$pen
  if (is.null(pen)) return("Penalty: unknown")
  choice <- if (length(x$pen_pen.choice) > 1) x$pen_pen.choice[2] else "custom"
  if (pen == "beta") {
    sprintf("Penalty: adaptive beta (choice %s), p = %.4g, q = %.4g",
            choice, x$p_q[1], x$p_q[2])
  } else if (pen == "norm") {
    sprintf("Penalty: adaptive normal (choice %s), mu = %.4g, sd = %.4g",
            choice, x$mu_std[1], x$mu_std[2])
  } else if (pen == "no") {
    "Penalty: none (pure L-moment fit)"
  } else {
    sprintf("Penalty: %s (fixed)", pen)
  }
}

convergence.label <- function(x) {
  if (is.null(x$convergence)) return("")
  if (identical(unname(x$convergence), 0)) "(converged)" else "(NOT converged)"
}

# Gumbel-scale Q-Q plot shared by the GEV11 fits --------------------------
qq.gum01 <- function(xdat, para, main, ...) {
  z <- trans.gum01(xdat, para = para, na.ok = TRUE)
  z <- sort(z[!is.na(z)])
  th <- -log(-log(ppoints(length(z))))
  plot(th, z, xlab = "Standard Gumbel quantiles",
       ylab = "Transformed observations", main = main, ...)
  abline(0, 1, lty = 2)
  invisible(NULL)
}

#' Methods for fitted GLmom objects
#'
#' @description
#' The estimation functions return classed objects: \code{\link{glme.gev}}
#' returns an object of class \code{"glme"}, \code{\link{glme.gev11}} of
#' class \code{"glme11"}, \code{\link{lme.gev11}} of class \code{"lme11"},
#' and \code{\link{ma.gev}} of class \code{"magev"}. The objects remain
#' ordinary lists, so all documented fields (e.g., \code{$para.glme}) stay
#' accessible as before.
#'
#' \code{print} shows the parameter estimates and the penalty settings;
#' \code{summary} additionally reports the auxiliary estimates stored in
#' the object; \code{plot} draws a quantile-quantile diagnostic: against
#' the fitted GEV distribution for \code{"glme"}, on the standard Gumbel
#' scale (after the transformation of Eq. 8 in the manuscript) for
#' \code{"glme11"} and \code{"lme11"}, and via \code{\link{magev.qqplot}}
#' for \code{"magev"}.
#'
#' @param x,object A fitted object.
#' @param ... Further arguments passed to plotting functions; ignored by
#'   \code{print} and \code{summary}.
#' @return \code{x} (or \code{object}), invisibly.
#'
#' @examples
#' data(haenam)
#' fit <- glme.gev(haenam$X1)
#' fit
#' summary(fit)
#' plot(fit)
#'
#' @name GLmom-methods
NULL

#' @rdname GLmom-methods
#' @export
print.glme <- function(x, ...) {
  cat("Generalized L-moment fit of the stationary GEV distribution\n")
  cat(penalty.label(x), "\n\n")
  cat("Estimates (para.glme):\n")
  print(round(x$para.glme, 4))
  cat("\nL-moment estimates (para.lme):\n")
  print(round(x$para.lme, 4))
  if (!is.null(x$nllh.glme))
    cat("\nPenalized criterion:", round(x$nllh.glme, 4),
        convergence.label(x), "\n")
  invisible(x)
}

#' @rdname GLmom-methods
#' @export
summary.glme <- function(object, ...) {
  print(object)
  if (!is.null(object$c0_c1_c2))
    cat("Hyperparameters (c0, c1, c2):",
        paste(object$c0_c1_c2, collapse = ", "), "\n")
  cat("Sample size:", length(object$data), "\n")
  invisible(object)
}

#' @rdname GLmom-methods
#' @export
plot.glme <- function(x, ...) {
  n <- length(x$data)
  th <- quagev(ppoints(n), vec2par(x$para.glme, "gev"))
  plot(th, sort(x$data), xlab = "Fitted GEV quantiles",
       ylab = "Empirical quantiles", main = "Q-Q plot of the GLME fit", ...)
  abline(0, 1, lty = 2)
  invisible(x)
}

#' @rdname GLmom-methods
#' @export
print.glme11 <- function(x, ...) {
  cat("Generalized L-moment fit of the nonstationary GEV11 model\n")
  cat("  mu(t) = mu0 + mu1 t,  sigma(t) = exp(sigma0 + sigma1 t)\n")
  cat(penalty.label(x), "\n\n")
  cat("Estimates (para.glme; mu1 and sigma1 held at the WLS pre-estimates):\n")
  print(round(x$para.glme, 4))
  if (!is.null(x$nllh.glme))
    cat("\nPenalized criterion:", round(x$nllh.glme, 4),
        convergence.label(x), "\n")
  if (!is.null(x$precis))
    cat("Precision of the L-moment equations:", signif(x$precis, 3), "\n")
  invisible(x)
}

#' @rdname GLmom-methods
#' @export
summary.glme11 <- function(object, ...) {
  print(object)
  cat("\nAuxiliary estimates:\n")
  aux <- rbind(WLS = object$para.wls, "WLS (original)" = object$strup.org)
  if (!is.null(object$para.gado)) aux <- rbind(aux, GN16 = object$para.gado)
  print(round(aux, 4))
  cat("\nStationary L-moment estimates (lme.sta):\n")
  print(round(object$lme.sta, 4))
  invisible(object)
}

#' @rdname GLmom-methods
#' @export
plot.glme11 <- function(x, ...) {
  qq.gum01(x$data, x$para.glme,
           main = "Gumbel-scale Q-Q plot of the GEV11 fit", ...)
  invisible(x)
}

#' @rdname GLmom-methods
#' @export
print.lme11 <- function(x, ...) {
  cat("L-moment fit of the nonstationary GEV11 model (Shin et al. 2025b)\n")
  cat("  mu(t) = mu0 + mu1 t,  sigma(t) = exp(sigma0 + sigma1 t)\n\n")
  cat("Estimates (lme.gev11):\n")
  print(round(x$lme.gev11, 4))
  cat("\nPrecision of the L-moment equations:", signif(x$precis, 3), "\n")
  invisible(x)
}

#' @rdname GLmom-methods
#' @export
summary.lme11 <- function(object, ...) {
  print(object)
  cat("Sample size:", length(object$data), "\n")
  invisible(object)
}

#' @rdname GLmom-methods
#' @export
plot.lme11 <- function(x, ...) {
  qq.gum01(x$data, x$lme.gev11,
           main = "Gumbel-scale Q-Q plot of the GEV11 fit", ...)
  invisible(x)
}

#' @rdname GLmom-methods
#' @export
print.magev <- function(x, ...) {
  cat("Model-averaged GEV quantile estimates\n")
  cat(sprintf("Weights: '%s' over %d candidate shape parameters\n\n",
              x$weight, x$run.numk))
  est <- rbind(MLE = x$qua.mle, LME = x$qua.lme, "Model average" = x$qua.ma)
  colnames(est) <- paste0("q", x$quant)
  print(round(est, 2))
  if (!is.null(x$fixw.se.ma)) {
    cat("\nSE of the model average (fixed weights):\n")
    print(round(x$fixw.se.ma, 2))
  }
  invisible(x)
}

#' @rdname GLmom-methods
#' @export
summary.magev <- function(object, ...) {
  print(object)
  cat("\nCandidate shape parameters (pick_xi):\n")
  print(round(object$pick_xi, 4))
  cat("\nWeights (w.ma):\n")
  print(round(object$w.ma, 4))
  invisible(object)
}

#' @rdname GLmom-methods
#' @export
plot.magev <- function(x, ...) {
  magev.qqplot(data = x$data, zx = x)
  invisible(x)
}
