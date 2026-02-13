# Shin et al. (2025b) compatibility functions
# Non-stationary GEV estimation using L-moments

#' Non-stationary GEV Parameter Estimation
#'
#' @description
#' Estimates parameters of a non-stationary Generalized Extreme Value (GEV)
#' distribution using the L-moment-based algorithm from Shin et al. (2025,
#' J. Korean Stat. Soc.).
#' This function combines L-moments, goodness-of-fit measures, and robust regression.
#'
#' This is a convenience wrapper around \code{\link{glme.gev11}} with \code{pen="no"},
#' providing compatibility with the original nsgev package interface.
#'
#' @param xdat A numeric vector of data to be fitted (e.g., annual maximum values).
#' @param ntry Number of attempts for optimization (default 20).
#' @param ftol Function tolerance for optimization (default 1e-6).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{para.prop} - The proposed L-moment based estimates
#'         (mu0, mu1, sigma0, sigma1, xi).
#'   \item \code{precis} - Precision of the optimization.
#' }
#'
#' @references
#' Shin, Y., Shin, Y. & Park, J.-S. (2025). Building nonstationary extreme value
#' model using L-moments. Journal of the Korean Statistical Society, 54, 947-970.
#' \doi{10.1007/s42952-025-00325-3}
#'
#' @seealso \code{\link{glme.gev11}} for the full GLME method with penalty functions,
#'   \code{\link{gado.prop_11}} for detailed estimation results.
#'
#' @author Yonggwan Shin, Seokkap Ko, Jihong Park, Yire Shin, Jeong-Soo Park
#'
#' @examples
#' data(Trehafod)
#' result <- nsgev(Trehafod$r1, ntry = 5)
#' print(result$para.prop)
#'
#' @export
nsgev <- function(xdat, ntry = 20, ftol = 1e-6) {
  z <- list()

  # Use glme.gev11 with no penalty (pure L-moment method)
  result <- glme.gev11(xdat, ntry = ntry, ftol = ftol,
                       init.rob = TRUE, pen = "no")

  # Return in nsgev format
  z$para.prop <- result$para.lme
  z$precis <- result$precis

  return(z)
}

#' Comprehensive Non-stationary GEV Estimation
#'
#' @description
#' Estimates parameters of a non-stationary GEV distribution using multiple methods:
#' Weighted Least Squares (WLS), GN16 method, and the proposed L-moment method
#' from Shin et al. (2025, J. Korean Stat. Soc.).
#'
#' This is a convenience wrapper around \code{\link{glme.gev11}} with \code{pen="no"},
#' providing compatibility with the original nsgev package interface.
#'
#' @param xdat A numeric vector of data to be fitted.
#' @param ntry Number of attempts for optimization (default 20).
#' @param ftol Function tolerance for optimization (default 1e-6).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{para.prop} - L-moment based estimates (proposed method).
#'   \item \code{para.gado} - GN16 method estimates.
#'   \item \code{para.wls} - Weighted least squares estimates.
#'   \item \code{strup.org} - Original non-stationary WLSE by Strup method.
#'   \item \code{lme.sta} - Stationary L-moment estimates.
#' }
#'
#' @references
#' Shin, Y., Shin, Y. & Park, J.-S. (2025). Building nonstationary extreme value
#' model using L-moments. Journal of the Korean Statistical Society, 54, 947-970.
#' \doi{10.1007/s42952-025-00325-3}
#'
#' @seealso \code{\link{glme.gev11}} for the full GLME method with penalty functions,
#'   \code{\link{nsgev}} for the simple interface.
#'
#' @author Yonggwan Shin, Seokkap Ko, Jihong Park, Yire Shin, Jeong-Soo Park
#'
#' @examples
#' data(Trehafod)
#' result <- gado.prop_11(Trehafod$r1, ntry = 5)
#' print(result$para.prop)
#' print(result$lme.sta)
#'
#' @export
gado.prop_11 <- function(xdat, ntry = 20, ftol = 1e-6) {
  z <- list()

  # Use glme.gev11 with no penalty (pure L-moment method)
  result <- glme.gev11(xdat, ntry = ntry, ftol = ftol,
                       init.rob = TRUE, pen = "no")

  # Return in gado.prop_11 format
  z$para.prop <- result$para.lme
  z$para.gado <- result$para.gado
  z$para.wls <- result$para.wls
  z$strup.org <- result$strup.org
  z$lme.sta <- result$lme.sta

  return(z)
}
