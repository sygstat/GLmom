# --------------------------------------------------------------
# Compute GN16 parameter estimates for the NS GEV11 model
# --------------------------------------------------------------

#' GN16 estimation for the non-stationary GEV11 model
#'
#' @description
#' Estimates the parameters of the non-stationary GEV11 model by the
#' quantile-based GN16 method (Gilleland and Katz, 2016; as adapted in
#' Shin et al., 2025b), with both the original and the modified
#' specification returned.
#'
#' @param xdat A numeric vector of data to be fitted.
#' @param rob If TRUE, uses robust regression in the internal steps
#'   (default FALSE).
#'
#' @return A list containing:
#' \itemize{
#'   \item para.gado.org - Original GN16 estimates
#'     (mu0, mu1, sigma0, sigma1, xi).
#'   \item para.gado.mdfy - Modified GN16 estimates.
#' }
#'
#' @references
#' Gilleland, E. & Katz, R. W. (2016). extRemes 2.0: An extreme value
#' analysis package in R. Journal of Statistical Software, 72(8), 1-39.
#'
#' Shin, Y., Shin, Y. & Park, J.-S. (2025). Building nonstationary extreme value
#' model using L-moments. Journal of the Korean Statistical Society, 54, 947-970.
#' \doi{10.1007/s42952-025-00325-3}
#'
#' @seealso \code{\link{glme.gev11}}, \code{\link{lme.gev11}},
#'   \code{\link{strup.gev11}}.
#'
#' @author Yonggwan Shin, Seokkap Ko, Jihong Park, Yire Shin, Jeong-Soo Park
#'
#' @examples
#' data(PhliuAgromet)
#' result <- GN16.gev11(PhliuAgromet$prec)
#' print(result$para.gado.org)
#'
#' @export
GN16.gev11 = function(xdat, rob=FALSE){

 year=seq(1,length(xdat))
 reg.dat=data.frame( cbind(year, xdat) )

 mu.init= lm(xdat~year, reg.dat)$coefficients

 orig.para=c(mu.init[1],mu.init[2], 1,.1,0)

 qlist= make.qmax.gev11(xdat, orig.para=orig.para, rob=rob)

 orig.para=c(mu.init[1],mu.init[2], qlist$sig0, qlist$sig1, 0)

 time.m.gev11(qmax=qlist$qmax, orig.para=orig.para,
              rob=rob)
}

#-----------------------------------------------------
#' Construct pseudo-maxima for the GN16 method
#'
#' @description Internal step of \code{\link{GN16.gev11}}: builds the
#' adjusted residual series (qmax) and estimates the log-scale trend
#' coefficients.
#'
#' @param xdat Numeric data vector.
#' @param orig.para Initial parameter vector (mu0, mu1, sigma0, sigma1, xi).
#' @param rob If TRUE, uses robust regression.
#' @return A list with sig0, sig1, and qmax.
#' @keywords internal
make.qmax.gev11 =function(xdat=NULL, orig.para=NULL, rob=FALSE)
{
  z=list()
  ns = length(xdat)
  year= seq(1,ns)

  res = xdat - (orig.para[1] + orig.para[2]* year)
  mres= mean(res)
  res.pr = abs(res - mres)

  lres.pr=log(res.pr)
  sig.dat=data.frame( cbind(year, lres.pr) )

  if(rob==FALSE){
    sig.lm= lm(lres.pr~year, sig.dat)$coefficients

  }else{
    sig.lm= lmrob(lres.pr~year, sig.dat)$coefficients
  }

  z$sig0 = sig.lm[1]
  z$sig1 = sig.lm[2]
  sigt = exp(z$sig0 + z$sig1* year)

  qmax= rep(NA, ns)

  for(i in 1:ns){
    if(z$sig1 >= 0){
      if(res[i] >= mres ) {
        qmax[i] = res[i] - sigt[i]
      }else{
        qmax[i] = res[i] + sigt[i]
      }
    }else{
      if(res[i] >= mres ) {
        qmax[i] = res[i] + sigt[i]
      }else{
        qmax[i] = res[i] - sigt[i]
      }
    }
  }
  z$qmax=qmax
  return(z)
}

#------------------------------------------------------------------
#' Time-varying moment step of the GN16 method
#'
#' @description Internal step of \code{\link{GN16.gev11}}: fits a
#' stationary GEV to the pseudo-maxima and converts the result to the
#' GEV11 parameterization (both original and modified specification).
#'
#' @param qmax Pseudo-maxima series from \code{\link{make.qmax.gev11}}.
#' @param orig.para Parameter vector (mu0, mu1, sigma0, sigma1, xi).
#' @param rob If TRUE, uses robust regression for the location trend.
#' @return A list with para.gado.org and para.gado.mdfy.
#' @keywords internal
time.m.gev11 = function(qmax=NULL, orig.para=NULL, rob=FALSE){

  z=list()
  ns=length(qmax)
  year=seq(1,ns)

  sig0=orig.para[3]
  sig1=orig.para[4]

  q.sta = pargev(lmoms(qmax, nmom=3), checklmom=FALSE)$para
  xi= q.sta[3]

  if( xi <= -0.5) xi = -0.499

  cd = sqrt( (xi^2) /( gamma(1+2*xi) - gamma(1+xi)^2 ) )
  alpha_t = exp(sig0 +sig1*year) * cd
  mu.gado = ( -( 1-gamma(1+ xi) )*alpha_t/xi
              + orig.para[1]+ orig.para[2]*year )

  #------------------------------
  nh=round((ns/2))
  mu.gado[nh-1] = mu.gado[nh-1] + 0.02
  mu.gado[nh+1] = mu.gado[nh+1] - 0.02

  mu.data= data.frame( cbind(year, mu.gado) )
  if(rob==TRUE){
    loc.gado =lmrob(mu.gado~year, mu.data)$coefficients

  }else{
    loc.gado =lm(mu.gado~year, mu.data)$coefficients
  }

  alpha0 = log(cd) + sig0

  z$para.gado.org= c(loc.gado, alpha0, sig1, xi)

  # modify GN16 ---------------

  alpha0_up = log(q.sta[2])- sig1*nh
  z$para.gado.mdfy = c(loc.gado, alpha0_up, sig1, xi)

  names(z$para.gado.org) <- names(z$para.gado.mdfy) <- c("mu0",
                                 "mu1","sigma0","sigma1","xi")
  return(z)
}
