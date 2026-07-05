# -----------------------------------------------------------------
# Building NS gev11 model using L-moments
# See the paper by Shin, Shin, and Park at JKSS (2025)
# -----------------------------------------------------------------

#' L-moment parameter estimation for the non-stationary GEV11 model
#'
#' @description
#' Estimates the parameters of the non-stationary GEV11 model
#' (mu(t) = mu0 + mu1*t, sigma(t) = exp(sigma0 + sigma1*t), constant xi)
#' by the pure L-moment method of Shin et al. (2025b, JKSS).
#' Starting from the WLS pre-estimate of \code{\link{strup.gev11}},
#' the three L-moment equations (matching the sample L-moments of the
#' Gumbel-transformed data with the population L-moments of the standard
#' Gumbel distribution) are solved by Broyden's method for
#' (mu0, sigma0, xi), and the best solution is chosen by a goodness-of-fit
#' measure.
#'
#' This function replaces \code{\link{gado.prop_11}} of GLmom v1.x
#' (renamed in v2.0.0, with slightly changed input/output).
#' It is also equivalent to \code{glme.gev11(xdat, pen="no")}.
#'
#' @param xdat A numeric vector of data to be fitted.
#' @param ntry Number of attempts (initial points) for the nonlinear solver
#'   (default 5).
#' @param ftol Function tolerance for the solver (default 1e-5).
#' @param show If TRUE, prints the precision and parameters of each try.
#'
#' @return A list containing:
#' \itemize{
#'   \item lme.gev11 - The L-moment estimates (mu0, mu1, sigma0, sigma1, xi).
#'     If no solution is found, the WLS pre-estimate is returned.
#'   \item precis - Precision (mean absolute residual of the L-moment
#'     equations) of the selected solution.
#'   \item data - The input data (used by the plot method).
#' }
#' The object has class \code{"lme11"}; see \code{\link{GLmom-methods}}.
#'
#' @references
#' Shin, Y., Shin, Y. & Park, J.-S. (2025). Building nonstationary extreme value
#' model using L-moments. Journal of the Korean Statistical Society, 54, 947-970.
#' \doi{10.1007/s42952-025-00325-3}
#'
#' @seealso \code{\link{glme.gev11}} for the penalized (GLME) version,
#'   \code{\link{strup.gev11}} for the WLS pre-estimation,
#'   \code{\link{GN16.gev11}} for the GN16 method,
#'   \code{\link{nsgev}} and \code{\link{gado.prop_11}} for the v1.x
#'   compatibility wrappers.
#'
#' @author Yonggwan Shin, Seokkap Ko, Jihong Park, Yire Shin, Jeong-Soo Park
#'
#' @examples
#' \donttest{
#' data(PhliuAgromet)
#' result <- lme.gev11(PhliuAgromet$prec, ntry = 5)
#' print(result$lme.gev11)
#' }
#'
#' @export
lme.gev11 = function(xdat, ntry=5, ftol=1e-5, show=FALSE){

  z <- list()
  model="gev11"
  init=matrix(NA, ntry,3)

  pretheta= strup.gev11(xdat, init.rob=TRUE)$strup.mdfy

  init[1,]= pretheta[c(1,3,5)]+0.02
  if(ntry >= 2){
   sinit = init.glme(xdat, ntry=ntry, pretheta=pretheta,
                     model=model)
   init[2:ntry,]= sinit[1:(ntry-1),]
  }

  precis=rep(1000, ntry)
  para.sel=matrix(NA,ntry,ncol=5)

  for(i in 1:ntry) {

    kt <- try(
            nleqslv(x=as.vector(init[i,1:3]),
                    fn= fun.lme.gev11,
                    method="Broyden",
                    control=list(ftol=ftol, xtol=ftol*10),
                    xdat=xdat, pretheta=pretheta),
            silent=TRUE)

    if(inherits(kt, "try-error")) next   # precis stays 1000

    precis[i]= mean(abs(kt$fvec))
    if(is.na(precis[i])) precis[i]= 1000

    if( precis[i] < ftol) {
      para.sel[i,1:5]= c(kt$x[1], pretheta[2],
                         kt$x[2], pretheta[4], kt$x[3])
    }

    if( abs( kt$termcd ) > 3 ) {
      precis[i]=1000
      para.sel[i,]=NA
    }

    if(show==TRUE){
     cat("--- itry, precision=",i, precis[i],"\n")
     cat("--- itry, para=",i, round(para.sel[i,1:5],4),"\n","\n")
    }
  } #end for

  sel = sel.para_all(xdat, para.sel, model,
                                     obj.fun=abs(precis))
  z$lme.gev11= sel$para
  z$precis = precis[sel$min.itry]

  if(z$precis > ftol) { z$lme.gev11 = pretheta
    message("no solution found for the L-moment method; ",
            "WLS pre-estimate returned") }

  names(z$lme.gev11) <- c("mu0","mu1","sigma0","sigma1","xi")
  z$data = xdat
  class(z) = "lme11"
  return(z)
}
