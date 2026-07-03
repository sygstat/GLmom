# --------------------------------------------------------------
# Random number generation and quantile functions for stationary
# and non-stationary GEV models
# --------------------------------------------------------------

#' Set parameters based on non-stationary model type
#'
#' @description Internal function that maps a flat parameter vector to named
#' location (mu), scale (sig), and shape (xi) components depending on the
#' model specification (GEV00, GEV10, GEV11, GEV20, GEV01).
#'
#' @param para Parameter vector whose length and meaning depend on \code{model}.
#' @param model Model type: "gev11", "gev10", "gev20", "gev01", or "gev00"/"gev".
#'
#' @return A list containing:
#' \describe{
#'   \item{mu}{Numeric vector c(mu0, mu1, mu2) for location}
#'   \item{sig}{Numeric vector c(sigma0, sigma1) for log-scale}
#'   \item{xi}{Shape parameter (scalar)}
#' }
#'
#' @keywords internal
set.para.model = function(para, model=NULL){

  z=list()
  mu0=para[1]
  mu1=para[2]

  if(model=='gev10'){
    sig0=para[3]
    xi=para[4]
    mu2=0
    sig1=0
  }else if(model=='gev20'){
    mu2=para[3]
    sig0=para[4]
    xi=para[5]
    sig1=0
  }else if(model=='gev11'){
    sig0=para[3]
    sig1=para[4]
    xi=para[5]
    mu2=0
  }else if(model=='gev00' | model=='gev'){
    mu1=mu2=sig1=0
    sig0=para[2]
    xi=para[3]
  }else if(model=='gev01'){
    mu1=mu2=0
    sig0=para[2]
    sig1=para[3]
    xi=para[4]
  }

  z$mu=c(mu0,mu1,mu2)
  z$sig=c(sig0,sig1)
  z$xi=xi
  return(z)
}

#---------------------------------------------------------
#' Random number generation for non-stationary GEV models
#'
#' @description
#' Generates a random sample from a stationary or non-stationary GEV model.
#' For non-stationary models the location and/or scale vary with time
#' t = 1, ..., nsample:
#' mu(t) = mu0 + mu1*t + mu2*t^2 and sigma(t) = exp(sigma0 + sigma1*t).
#'
#' @param nsample Sample size (number of time points).
#' @param para Parameter vector; its length and meaning depend on
#'   \code{model} (see \code{\link{quagev.NS}}). For "gev11":
#'   (mu0, mu1, sigma0, sigma1, xi).
#' @param model Model type: "gev", "gev00", "gev01", "gev10", "gev11",
#'   or "gev20".
#'
#' @return A numeric vector of length \code{nsample}.
#'
#' @references
#' Shin, Y., Shin, Y. & Park, J.-S. (2025). Building nonstationary extreme value
#' model using L-moments. Journal of the Korean Statistical Society, 54, 947-970.
#' \doi{10.1007/s42952-025-00325-3}
#'
#' @seealso \code{\link{quagev.NS}} for the corresponding quantile function,
#'   \code{\link{glme.gev11}}, \code{\link{lme.gev11}}.
#'
#' @author Yonggwan Shin, Seokkap Ko, Jihong Park, Yire Shin, Jeong-Soo Park
#'
#' @examples
#' set.seed(1)
#' # Non-stationary GEV11 sample of size 50
#' x <- ran.gev_all(50, para = c(100, 0.5, 3, 0.005, -0.2), model = "gev11")
#' plot(x, type = "l")
#'
#' @export
ran.gev_all= function(nsample, para, model=NULL){

  rns=rep(NA, nsample)

  # model=c("gev","gev01","gev10", "gev11", "gev20")

  sp=set.para.model(para,model)

  for (t in 1:nsample){
    mut=sp$mu[1]+ sp$mu[2]*t +sp$mu[3]*(t^2)
    sigt=exp(sp$sig[1] +sp$sig[2]*t)

    savens=vec2par(c(mut,sigt,sp$xi),'gev')
    rns[t]=rlmomco(1,savens)
  }

  return(rns)
}

#----------------------------------------------------
#' Quantile function for non-stationary GEV models
#'
#' @description
#' Calculates quantiles for non-stationary GEV models including GEV11,
#' GEV10, GEV20, and stationary GEV00.
#'
#' @param f Probability (or vector of probabilities) for quantile calculation.
#' @param para Parameter vector. For GEV11: (mu0, mu1, sigma0, sigma1, xi).
#' @param nsample Number of time points (sample size).
#' @param model Model type: "gev11", "gev10", "gev20", "gev01", or "gev00"/"gev".
#'
#' @return A matrix of quantiles (nsample x length(f)) or a vector if f is scalar.
#'
#' @references
#' Shin, Y., Shin, Y., Park, J. & Park, J.-S. (2025). Generalized method of
#' L-moment estimation for stationary and nonstationary extreme value models.
#' arXiv preprint arXiv:2512.20385. \doi{10.48550/arXiv.2512.20385}
#'
#' @seealso \code{\link{glme.gev11}} for non-stationary GEV estimation,
#'   \code{\link{glme.gev}} for stationary GEV estimation,
#'   \code{\link{ran.gev_all}} for random number generation.
#'
#' @author Yonggwan Shin, Seokkap Ko, Jihong Park, Yire Shin, Jeong-Soo Park
#'
#' @examples
#' # GEV11 model: time-varying quantiles
#' para <- c(84.55, 1.03, 2.91, 0.009, -0.08)  # mu0, mu1, sigma0, sigma1, xi
#' q99 <- quagev.NS(f = 0.99, para = para, nsample = 53, model = "gev11")
#' print(q99)
#'
#' @export
quagev.NS= function(f=NULL, para=NULL, nsample=NULL, model=NULL){

  if(model=='gev00' | model=='gev'){
    zpT = quagev(f,vec2par(para,'gev'))
    return(zpT)

  }else{  # if(model !='gev00')

    ns=nsample
    year=seq(1:ns)
    numq=length(f)
    zpT= matrix(NA,ns,numq)

    sp=set.para.model(para,model)

    xi= sp$xi
    zpc= (1- (-log(f))^xi) /xi
    vec = sp$mu[1] + sp$mu[2]*year[1:ns] +sp$mu[3]*(year[1:ns]^2)

    for(iq in 1:numq){
      zpT[1:ns,iq] = vec[1:ns]+ zpc[iq]* exp( sp$sig[1]
                                               +sp$sig[2]*year[1:ns] )
    }
    if(numq==1) zpT= as.vector(zpT)
    zpT
  }
}

#----------------------------------------------------
#' Quantile function for GEV11 model
#'
#' @description Internal wrapper that computes time-varying quantiles for the
#' GEV11 model given a return period.
#'
#' @param Tp Return period.
#' @param para Parameter vector (mu0, mu1, sigma0, sigma1, xi).
#' @param year Time vector.
#' @return Numeric vector of quantile values (one per time point).
#' @keywords internal
qns.gev11= function(Tp=NULL, para=NULL, year=NULL){

  model="gev11"
  f=1-1/Tp
  quagev.NS(f,para,nsample=length(year),model)
}

#----------------------------------------------------
#' Quantile function for non-stationary GEV (return period input)
#'
#' @description Internal function that computes time-varying quantiles for
#' non-stationary GEV models given a return period.
#'
#' @param Tp Return period (scalar).
#' @param para Parameter vector for the non-stationary model.
#' @param year Numeric vector of time points.
#' @param model Model type: "gev11", "gev10", "gev20", or "gev01".
#'
#' @return Numeric vector of quantile values (one per time point).
#'
#' @keywords internal
qns.gev_all= function(Tp=NULL, para=NULL, year=NULL, model=NULL){

  nsample=length(year)
  ns=nsample
  zpT=rep(NA, nsample)
  year2=year^2

  sp=set.para.model(para,model)

  xi= sp$xi
  zpc= (1- ( -log(1-(1/Tp) ) )^xi ) /xi
  zpT[1:ns] = sp$mu[1] + sp$mu[2]*year[1:ns] +sp$mu[3]*year2[1:ns]
  zpT[1:ns] = zpT[1:ns]+ zpc* exp( sp$sig[1] +sp$sig[2]*year[1:ns] )

  return(zpT)
}
