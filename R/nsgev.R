# Shin et al. (2025b) compatibility functions
# Non-stationary GEV estimation using L-moments

#' Non-stationary GEV Parameter Estimation
#'
#' @description
#' Estimates parameters of a non-stationary Generalized Extreme Value (GEV)
#' distribution using the L-moment-based algorithm from Shin et al. (2025,
#' J. Korean Stat. Soc.).
#'
#' This is a convenience wrapper around \code{\link{lme.gev11}}, providing
#' compatibility with the original nsgev package interface (and with
#' GLmom v1.x).
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
#' @seealso \code{\link{lme.gev11}} for the recommended interface,
#'   \code{\link{glme.gev11}} for the full GLME method with penalty functions,
#'   \code{\link{gado.prop_11}} for detailed estimation results.
#'
#' @author Yonggwan Shin, Seokkap Ko, Jihong Park, Yire Shin, Jeong-Soo Park
#'
#' @examples
#' data(PhliuAgromet)
#' result <- nsgev(PhliuAgromet$prec, ntry = 5)
#' print(result$para.prop)
#'
#' @export
nsgev <- function(xdat, ntry = 20, ftol = 1e-6) {
  z <- list()

  result <- lme.gev11(xdat, ntry = ntry, ftol = ftol)

  # Return in nsgev format
  z$para.prop <- result$lme.gev11
  z$precis <- result$precis

  return(z)
}

#' Comprehensive Non-stationary GEV Estimation (deprecated)
#'
#' @description
#' Estimates parameters of a non-stationary GEV distribution using multiple
#' methods: Weighted Least Squares (WLS), GN16 method, and the proposed
#' L-moment method from Shin et al. (2025, J. Korean Stat. Soc.).
#'
#' As of GLmom v2.0.0 this function is deprecated: it was renamed to
#' \code{\link{lme.gev11}} (for the proposed method), and the auxiliary
#' methods are available as \code{\link{strup.gev11}} and
#' \code{\link{GN16.gev11}}. This wrapper reassembles the v1.x output
#' format from those functions and will be removed in a future release.
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
#' @seealso \code{\link{lme.gev11}} for the recommended interface,
#'   \code{\link{strup.gev11}}, \code{\link{GN16.gev11}},
#'   \code{\link{glme.gev11}} for the full GLME method,
#'   \code{\link{nsgev}} for the simple interface.
#'
#' @author Yonggwan Shin, Seokkap Ko, Jihong Park, Yire Shin, Jeong-Soo Park
#'
#' @examples
#' \donttest{
#' data(PhliuAgromet)
#' # deprecated wrapper; suppressWarnings() silences the deprecation notice
#' result <- suppressWarnings(gado.prop_11(PhliuAgromet$prec, ntry = 5))
#' print(result$para.prop)
#' print(result$lme.sta)
#' }
#'
#' @export
gado.prop_11 <- function(xdat, ntry = 20, ftol = 1e-6) {
  .Deprecated("lme.gev11")
  z <- list()

  result <- lme.gev11(xdat, ntry = ntry, ftol = ftol)
  z$para.prop <- result$lme.gev11

  strup <- strup.gev11(xdat, init.rob = TRUE)
  z$para.wls <- strup$strup.mdfy
  z$strup.org <- strup$strup.para

  gado <- GN16.gev11(xdat)
  z$para.gado <- gado$para.gado.org

  z$lme.sta <- pargev(lmoms(xdat, nmom = 3))$para
  names(z$lme.sta) <- c("mu_sta", "sigma_sta", "xi_sta")

  return(z)
}
