# ------------------------------------------------------------------
# Estimate parameters of NS GEV11 by the Weighted Least Squares
# method of Strupczewski and Kaczmarek (2001)
# ------------------------------------------------------------------

#' Weighted least squares estimation for the non-stationary GEV11 model
#'
#' @description
#' Estimates the parameters of the non-stationary GEV11 model by the
#' weighted least squares (WLS) method of Strupczewski and Kaczmarek (2001),
#' with the final-specification modification described in Shin et al.
#' (2025b, JKSS). The result is used as the pre-estimate (starting value)
#' by \code{\link{glme.gev11}} and \code{\link{lme.gev11}}.
#'
#' @param xdat A numeric vector of data to be fitted.
#' @param init.rob If TRUE (default), the initial location trend is fitted
#'   by robust regression (\code{\link[robustbase]{lmrob}}); otherwise by
#'   ordinary least squares.
#' @param wls.rob If TRUE, the log-scale regression uses robust regression
#'   (default FALSE).
#' @param strup.fin.rob If TRUE (default), the final location trend
#'   re-specification uses robust regression.
#'
#' @return A list containing:
#' \itemize{
#'   \item strup.para - Original WLS estimates (mu0, mu1, sigma0, sigma1, xi).
#'   \item strup.mdfy - Modified (re-specified) WLS estimates.
#' }
#'
#' @references
#' Strupczewski, W. G. & Kaczmarek, Z. (2001). Non-stationary approach to
#' at-site flood frequency modelling II. Weighted least squares estimation.
#' Journal of Hydrology, 248(1-4), 143-151.
#'
#' Shin, Y., Shin, Y. & Park, J.-S. (2025). Building nonstationary extreme value
#' model using L-moments. Journal of the Korean Statistical Society, 54, 947-970.
#' \doi{10.1007/s42952-025-00325-3}
#'
#' @seealso \code{\link{glme.gev11}}, \code{\link{lme.gev11}},
#'   \code{\link{GN16.gev11}}.
#'
#' @author Yonggwan Shin, Seokkap Ko, Jihong Park, Yire Shin, Jeong-Soo Park
#'
#' @examples
#' data(PhliuAgromet)
#' result <- strup.gev11(PhliuAgromet$prec)
#' print(result$strup.mdfy)
#'
#' @export
strup.gev11 = function(xdat, init.rob=TRUE,
                       wls.rob=FALSE, strup.fin.rob=TRUE){
  w=list()
  ns=length(xdat)
  year=seq(1,ns)

  reg.dat=data.frame( cbind(year, xdat) )

  if(init.rob==FALSE){
    mu.init= lm(xdat~year, reg.dat)$coefficients  # regression
  }else{
    mu.init= lmrob(xdat~year, reg.dat)$coefficients  # robust regression
  }
  res= xdat -(mu.init[1]+mu.init[2]*year)

  wls= wls.gev11(xdat, res=res, rob=wls.rob)  # log-scale regression

  wls$res = res/exp(wls$sig[1] + wls$sig[2]*year)
  strup.sta = pargev(lmoms(wls$res), checklmom=FALSE)$para

  strup.para = c(mu.init[1:2], wls$sig, strup.sta[3])
  w$strup.para = strup.para

  #------to specify final parameter values ----------
  # See JKSS(2025) paper by Shin, Shin, Park

  yt=rep(0,ns+1)
  mu_st=strup.sta[1]
  year2=seq(0,ns)

  for (ka in 1:(ns+1) ) {
    yt[ka]= ( strup.para[1] + strup.para[2]*(ka-1)
            + mu_st*exp(strup.para[3]+ strup.para[4]*(ka-1)) )
  }
  nh=round((ns/2))
  yt[nh-1] = yt[nh-1] + 0.02;   yt[nh-2] = yt[nh-2] - 0.02
  yt[nh+1] = yt[nh+1] - 0.01;   yt[nh+2] = yt[nh+2] + 0.01

  reg.dat=data.frame( cbind(year2, yt) )

  if(strup.fin.rob==FALSE){
    mu.f= lm(yt ~ year2, reg.dat)$coefficients
  }else{
    mu.f= lmrob(yt ~ year2, reg.dat)$coefficients
  }
  sigmaf_0 = strup.para[3] + log(strup.sta[2])

  w$strup.mdfy= c(mu.f, sigmaf_0, strup.para[4], strup.sta[3])
  names(w$strup.para) <- names(w$strup.mdfy) <- c("mu0",
                               "mu1","sigma0","sigma1","xi")
  w
}

#------------------------------------------------------------------
#' Log-scale regression step of the Strupczewski WLS method
#'
#' @description Internal step of \code{\link{strup.gev11}}: regresses the
#' log absolute residuals on time to estimate the log-scale coefficients
#' (sigma0, sigma1).
#'
#' @param xdat Numeric data vector.
#' @param res Residuals from the location trend fit.
#' @param rob If TRUE, uses robust regression (setting "KS2014").
#' @return A list with element \code{sig} (the two coefficients).
#' @keywords internal
wls.gev11= function(xdat, res=NULL, rob=NULL){

  z=list()
  ns=length(res)
  year=seq(1,ns)

  lres.pr=log(abs(res))
  sig.dat=data.frame( cbind(year, lres.pr) )

  if(rob==FALSE){
    z$sig= lm(lres.pr~year, sig.dat)$coefficients
  }else{
    z$sig= lmrob(lres.pr~year, data=sig.dat,
                 setting="KS2014")$coefficients
  }
  return(z)
}
