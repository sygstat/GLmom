# R code for Non-Stationary GEV11 model using generalized L-moments
# GEV11: mu(t) = mu0 + mu1*t, sigma(t) = exp(sigma0 + sigma1*t), xi constant
# v2.0.0: rewritten following the revised code of Shin et al. (Apr 2026)

#----------------------------------------------------------------------
#' Generalized L-moments estimation for non-stationary GEV11 model
#'
#' @description
#' This function estimates parameters of the non-stationary GEV11 model
#' where mu(t) = mu0 + mu1*t and sigma(t) = exp(sigma0 + sigma1*t),
#' by the generalized method of L-moment estimation (GLME).
#'
#' Starting from the WLS pre-estimate of \code{\link{strup.gev11}},
#' the penalized objective (generalized L-moment distance based on the
#' asymptotic covariance of the sample L-moments, plus penalty) is
#' minimized over (mu0, sigma0, xi) with mu1 and sigma1 held at their
#' pre-estimated values. With \code{pen="no"} the result equals the pure
#' L-moment method \code{\link{lme.gev11}} of Shin et al. (2025b).
#'
#' @param xdat A numeric vector of data to be fitted.
#' @param ntry Number of attempts for parameter estimation (default 5).
#' @param ftol Tolerance for convergence (default 1e-6).
#' @param opt.choose Selection criterion among the multi-start solutions:
#'   "nllh" (default, penalized negative log-likelihood) or "gof"
#'   (goodness-of-fit by the number-of-exceedances measure).
#' @param pen Type of penalty function: "norm", "beta" (default), "ms",
#'   "park", "cannon", "cd", or "no".
#' @param pen.choice Choice number for penalty hyperparameters (default 1).
#'   See \code{\link{glme.gev}} for the preset values. Set to NULL to use
#'   \code{p, c1, c2} (beta) or \code{mu, std} (norm) directly.
#' @param mu Mean for normal penalty (default -0.55).
#' @param std Std for normal penalty (default 0.2).
#' @param p Shape for beta penalty (default 6).
#' @param c1 Scaling for beta penalty (default 5).
#' @param c2 Limit for beta penalty (default 2).
#' @param c0 Half-width of the adaptive beta penalty support (default 0.35).
#' @param q Optional fixed second shape parameter for the beta penalty.
#' @param show If TRUE, prints the objective value of each try.
#' @param init.rob Use robust regression for the WLS initialization
#'   (default TRUE).
#' @param glme.pre Deprecated (ignored). The WLS pre-estimate is always used.
#'
#' @return A list containing:
#' \itemize{
#'   \item para.glme - Proposed GLME estimates (mu0, mu1, sigma0, sigma1, xi).
#'   \item nllh.glme - Penalized negative log-likelihood of the solution.
#'   \item convergence - 0 if converged, 5 otherwise.
#'   \item pen, pen_pen.choice - Penalty method (and choice) used.
#'   \item p_q, c1_c2, c0_c1_c2 - (for beta) hyperparameters used.
#'   \item mu_std - (for norm) hyperparameters used.
#'   \item strup.org - WLSE by the Strupczewski method (\code{strup.para}).
#'   \item para.wls - Modified WLS estimates (\code{strup.mdfy}).
#'   \item para.gado - GN16 method estimates.
#'   \item lme.sta - Stationary L-moment estimates.
#'   \item para.lme - Pure L-moment estimates; returned only when
#'     \code{pen="no"} (then identical to para.glme). For other penalties
#'     use \code{\link{lme.gev11}} directly.
#'   \item precis - (only when \code{pen="no"}) precision of the L-moment
#'     equation solution.
#'   \item data - The input data (used by the plot method).
#' }
#' The object has class \code{"glme11"}; see \code{\link{GLmom-methods}}.
#'
#' @references
#' Shin, Y., Shin, Y., Park, J. & Park, J.-S. (2025). Generalized method of
#' L-moment estimation for stationary and nonstationary extreme value models.
#' arXiv preprint arXiv:2512.20385. \doi{10.48550/arXiv.2512.20385}
#'
#' Shin, Y., Shin, Y. & Park, J.-S. (2025). Building nonstationary extreme value
#' model using L-moments. Journal of the Korean Statistical Society, 54, 947-970.
#' \doi{10.1007/s42952-025-00325-3}
#'
#' @seealso \code{\link{glme.gev}} for stationary GEV estimation,
#'   \code{\link{lme.gev11}} for the pure L-moment method (no penalty),
#'   \code{\link{strup.gev11}} for the WLS pre-estimation,
#'   \code{\link{GN16.gev11}} for the GN16 method,
#'   \code{\link{quagev.NS}} for non-stationary quantile computation.
#'
#' @author Yonggwan Shin, Seokkap Ko, Jihong Park, Yire Shin, Jeong-Soo Park
#'
#' @examples
#' data(PhliuAgromet)
#' x <- PhliuAgromet$prec
#'
#' \donttest{
#' # Estimate non-stationary GEV11 parameters
#' result <- glme.gev11(x, ntry = 5)
#' print(result$para.glme)  # Proposed GLME estimates
#' }
#'
#' @export
glme.gev11= function(xdat, ntry=5, ftol=1e-6,
                     opt.choose="nllh", pen='beta', pen.choice=1,
                     mu= -0.55, std=0.2, p=6, c1=5, c2=2,
                     c0=0.35, q=NULL, show=FALSE,
                     init.rob=TRUE, glme.pre=NULL){

  if(!is.null(glme.pre)){
    warning("'glme.pre' is deprecated and ignored; ",
            "the WLS (Strupczewski) pre-estimate is always used.")
  }

  z <- list()
  model='gev11'
  name_gev11_ns  =c("mu0","mu1","sigma0","sigma1","xi")
  name_gev00_sta =c("mu_sta","sigma_sta","xi_sta")

  ch= check.penalty(pen=pen, pen.choice=pen.choice, p=p,
                   c1=c1,c2=c2,mu=mu,std=std)
  pen= ch$pen
  if(!is.null(pen.choice)){
   p= ch$p ; c1 =ch$c1;  c2= ch$c2
   mu= ch$mu ; std= ch$std
  }

#---- proposed method with pen="no" or "norm" or "beta" or etc.--
# z$para.glme = proposed glme

  strup = strup.gev11(xdat, init.rob=init.rob)
  pretheta= strup$strup.mdfy

  z = opt.glme.gev11(xdat, ntry=ntry, ftol=ftol,
                    pretheta=pretheta, model=model,
                    pen=pen, mu=mu, std=std,
                    p=p, c0=c0, c1=c1, c2=c2, q=q,
                    opt.choose=opt.choose, show=show)

  if(z$convergence > 0) { z$para.glme = pretheta
     message("no result was found for proposed glme; ",
             "WLS pre-estimate returned") }

  names(z$para.glme) <- name_gev11_ns
  z$pen = pen                          # kept for v1.x compatibility
  z$pen_pen.choice= c(pen, pen.choice)

  if(pen=='beta'){
    ww= pk.beta(para=z$para.glme[c(1,3,5)],
                lme.center=pretheta[c(1,3,5)],
                p=p, c0=c0, c1=c1, c2=c2, q=q)
    z$p_q= c(ww$p,ww$q)
    z$c1_c2= c(c1,c2)                  # kept for v1.x compatibility
    z$c0_c1_c2=c(c0,c1,c2)

  }else if(pen=="norm"){
    z$mu_std =c(mu,std)
  }

  #---- v1.x compatibility fields (inexpensive) ---------------------
  z$strup.org = strup$strup.para       # WLSE by Strupczewski
  z$para.wls  = strup$strup.mdfy       # modified WLSE
  names(z$strup.org) <- names(z$para.wls) <- name_gev11_ns

  gado <- try(GN16.gev11(xdat), silent=TRUE)
  if(!inherits(gado, "try-error")){
    z$para.gado = gado$para.gado.org   # GN16 original estimates
    names(z$para.gado) <- name_gev11_ns
  }

  z$lme.sta = pargev(lmoms(xdat,nmom=3))$para   # stationary L-ME
  names(z$lme.sta) <- name_gev00_sta

  if(pen=="no") z$para.lme = z$para.glme

  z$data = xdat
  class(z) = "glme11"
  return(z)
}

#---------------------------------------------------------------------
#' GLME objective function for GEV11 model (mu0, sigma0, xi optimization)
#'
#' @param a Parameter vector (mu0, sigma0, xi).
#' @param xdat Data vector.
#' @param newtheta Full parameter vector (mu0, mu1, sigma0, sigma1, xi).
#' @param covinv Inverse covariance matrix of the sample L-moments.
#' @param lcovdet Log determinant of the covariance matrix.
#' @param pen Penalty type.
#' @param mu Normal penalty mean.
#' @param std Normal penalty std.
#' @param p Beta penalty shape.
#' @param c0 Beta penalty support half-width.
#' @param c1 Beta penalty scaling.
#' @param c2 Beta penalty limit.
#' @param q Optional fixed second beta shape parameter.
#' @return Penalized negative log-likelihood value.
#' @author Yonggwan Shin, Seokkap Ko, Jihong Park, Yire Shin, Jeong-Soo Park
#' @keywords internal
nllh.glme.gev11 = function(a, xdat=xdat, newtheta=newtheta,
                           covinv=covinv, lcovdet=lcovdet,
                           pen=pen, mu=mu, std=std,
                           p=p, c0=c0, c1=c1, c2=c2, q=q)
{
  # updated 18Apr26
  ldist= fun.lme.gev11(a, xdat=xdat, pretheta=newtheta)

  if(any(ldist >= 10^5)) return(1e6)

  gld= t(ldist) %*% covinv %*% ldist

  nllh.norm =  gld/2 + (3/2)*log(2*pi) + lcovdet

  pk.ns= penalty.fun(par=a, mu=mu, std=std, lme=newtheta[c(1,3,5)],
                     pen=pen, p=p, c0=c0, c1=c1, c2=c2, q=q)
  as.numeric(nllh.norm + pk.ns)
}

#-----------------------------------------------------------------
#' L-moment distance function for GEV11 model
#'
#' @description Computes the vector of differences between the population
#' L-moments of the standard Gumbel distribution and the sample L-moments
#' of the data transformed by \code{\link{trans.gum01}} under the current
#' GEV11 parameters. Solving this system equal to zero gives the L-moment
#' estimates (Shin et al., 2025b).
#'
#' @param a Numeric vector (mu0, sigma0, xi) to optimize.
#' @param xdat Numeric vector of data.
#' @param pretheta Full parameter vector (mu0, mu1, sigma0, sigma1, xi)
#'   from pre-estimation; mu1 and sigma1 are held fixed.
#'
#' @return Numeric vector of length 3 of L-moment equation residuals.
#' @keywords internal
fun.lme.gev11 <- function(a, xdat=xdat, pretheta=pretheta)
{
  if( abs(a[3]) >= 1) return(rep(10^6,3))
  par=c(a[1],pretheta[2],a[2],pretheta[4],a[3])

  newg2= trans.gum01(xdat,para=par, na.ok=TRUE)
  newg=newg2[!is.na(newg2)]

  if( length(newg) < length(xdat)/2 ) return(rep(10^6,3))

  # population L-moments of the standard Gumbel:
  # lambda1 = Euler's constant, lambda2 = log(2), lambda3 = 0.169925*log(2)

  lgum= lmoms.md.park(newg, nmom=3, mtrim=FALSE, no.stop=TRUE)
  if(lgum$ifail == 1) return(rep(10^6,3))

  c(0.57721566, 0.69314718, 0.11778304) - lgum$lambdas[1:3]
}

#-----------------------------------------------------------------------
#' Multi-start penalized optimization for the GEV11 model
#'
#' @description Internal driver of \code{\link{glme.gev11}}. For
#' \code{pen="no"} it delegates to \code{\link{lme.gev11}}. Otherwise it
#' minimizes \code{\link{nllh.glme.gev11}} from several starting points,
#' using the fixed asymptotic covariance matrix of the first three sample
#' L-moments (computed for n=300 and rescaled by 300/n; see Eq. (26) of
#' Shin et al., 2025).
#'
#' @param xdat Numeric data vector.
#' @param ntry Number of starting points.
#' @param ftol Convergence tolerance.
#' @param show If TRUE, prints progress.
#' @param pretheta Pre-estimated parameter vector.
#' @param model Model name ("gev11").
#' @param pen,mu,std,p,q,c0,c1,c2 Penalty settings.
#' @param opt.choose "nllh" or "gof".
#' @return A list with para.glme, nllh.glme, convergence
#'   (and precis when \code{pen="no"}).
#' @keywords internal
opt.glme.gev11= function(xdat, ntry=5, ftol=1e-5, show=FALSE,
                        pretheta=NULL, model="gev11",
                        pen='beta', mu=NULL, std=NULL, p=2, q=NULL,
                        c0=0.35, c1=3, c2=2, opt.choose="nllh"){
  zm=list()

 if(pen=="no"){

   lme= lme.gev11(xdat, ntry=ntry, ftol=ftol)
   zm$precis = lme$precis

   if(lme$precis <= ftol){
     zm$para.glme = lme$lme.gev11
     zm$convergence =0
   }else{
     zm$para.glme = pretheta
     zm$convergence =5
   }
   return(zm)
 }

 # pen != "no" : perform glme -----------------------------------

  cov_300 = matrix(c( 0.0053986209, 0.0015941445, 0.0004516406,
                      0.0015941445, 0.0012886186, 0.0005222434,
                      0.0004516406, 0.0005222434, 0.0005117520),
                   3,3)

  cov_n0 = cov_300 *300/length(xdat)

  covinv= solve(cov_n0)
  lcovdet=log(det(cov_n0))

  newtheta = pretheta
  init = matrix(0, nrow=ntry, ncol=3)
  init[1,] = pretheta[c(1,3,5)]

  if(ntry >= 2){
    max_pk= find_max_beta.pk(newtheta[c(1,3,5)],p=p,
                             c0=c0,c1=c1,c2=c2)
    init[2,]= c(newtheta[c(1,3)]+0.2, max_pk)
  }
  if(ntry >= 3){
    sinit = init.glme(xdat, ntry=ntry, model="gev11",
                            pretheta=pretheta)
    init[3:ntry,]= sinit[1:(ntry-2),]
  }

  conv= my.nllh= rep(1e6,ntry)
  para.sel=matrix(NA,ntry,ncol=5)

  for(i in 1:ntry){

    kt <- try(
          optim(par=as.vector(init[i,1:3]),
                fn= nllh.glme.gev11,
                method="BFGS",
                control=list(abstol=ftol, reltol=ftol,
                               maxit=80),
                xdat=xdat, newtheta=newtheta, covinv=covinv,
                lcovdet=lcovdet, pen=pen, mu=mu,
                std=std, p=p, c0=c0, c1=c1, c2=c2, q=q),
      silent=TRUE)

    if(inherits(kt, "try-error")){
      my.nllh[i]= 10^6
      conv[i]= 5
      next
    }

    if(show==TRUE) cat("---itry, fvec=",i,kt$value,"\n")

    if( kt$convergence != 0) {
      my.nllh[i]=10^6
      conv[i] =5
    }else{
      my.nllh[i] = kt$value
      conv[i]= 0
      para.sel[i,1:5]=c(kt$par[1], newtheta[2],
                        kt$par[2], newtheta[4], kt$par[3])
    } #if
  } #for

  if(all(my.nllh==10^6)) {
    message("-- No solution found in optim for glme.gev11")
    zm$para.glme = pretheta
    zm$convergence= 5
    return(zm)
   }

  if(opt.choose=="nllh"){

      id= which.min(my.nllh)
      zm$para.glme = para.sel[id,]
      zm$nllh.glme = my.nllh[id]
      zm$convergence = conv[id]

  }else if(opt.choose=="gof"){

      sel = sel.para_all(xdat,    #para est. SSP JKSS(2025)
                         para.sel, model, obj.fun=my.nllh)

      zm$para.glme = sel$para
      zm$nllh.glme = my.nllh[sel$min.itry]
      zm$convergence = conv[sel$min.itry]
  } #if opt

 return(zm)
}

#-------------------------------------------------------------
#' Find the shape parameter maximizing the beta penalty
#'
#' @description Finds the value of xi (over a grid within the adaptive
#' beta penalty support) which maximizes the beta penalty function.
#' Used as one of the initial points for optimization.
#'
#' @param lme.center Center vector (mu0, sigma0, xi) of the penalty.
#' @param p,c0,c1,c2 Beta penalty hyperparameters.
#' @return The xi value maximizing the penalty.
#' @keywords internal
find_max_beta.pk=function(lme.center,p=p,
                          c0=c0,c1=c1,c2=c2){

  lower= max(-1, lme.center[3]-c0)
  upper= min(0.3, lme.center[3]+c0)
  xpk= seq(lower,upper,by=(upper-lower)/70)
  pk=rep(NA,length(xpk))

  for(i in 1:length(xpk)){
    a= c(0,1,xpk[i])
    pk[i]= pk.beta(para= a, lme.center=lme.center,
                   p=p, c0=c0,c1=c1,c2=c2)$pk.one
  }
  xpk[which.max(pk)]
}

#-----------------------------------------------------------
#' Select parameters using goodness-of-fit criterion
#'
#' @description Internal function that selects the best parameter set among
#' candidate solutions using a goodness-of-fit measure based on the expected
#' vs. observed number of exceedances (see Shin et al., 2025b, JKSS).
#'
#' @param xdat Numeric vector of data.
#' @param para.sel Matrix of candidate parameter vectors (one per row).
#' @param model Model name string (e.g., "gev11").
#' @param obj.fun Numeric vector of objective values used for tie-breaking.
#'
#' @return A list containing:
#' \describe{
#'   \item{para}{Selected parameter vector}
#'   \item{min.itry}{Index of selected try}
#'   \item{gof}{Goodness-of-fit values}
#'   \item{obj.fun}{Objective function values}
#' }
#' @keywords internal
sel.para_all =function(xdat, para.sel=NULL, model=NULL,
                       obj.fun=NULL){
  z=list()
  upara.sel = para.sel

  gof=rep(NA,nrow(upara.sel))
  ns=length(xdat)
  npar=ncol(upara.sel)

  vecT=c(5,10,20,40,60)
  if(ns >= 100) vecT=c(5,10,20,40,80,120)
  if(ns <= 30) vecT=c(5,10,20,40)

  for(i in 1:nrow(upara.sel) ){
    par.vec = as.vector(upara.sel[i,1:npar])
    gof[i]  = gof.ene_all(xdat, vecT, par.vec, model)
  }

  if(length(unique(gof))==1) {
    z$para =upara.sel[which.min(obj.fun),]
    z$min.itry = which.min(obj.fun)
  }else{
    z$para =upara.sel[which.min(gof),]
    z$min.itry = which.min(gof)
  }

  z$gof=gof
  z$obj.fun= obj.fun
  return(z)
}

#-----------------------------------------------------------
#' Goodness-of-fit by number of exceedances
#'
#' @description Internal function computing a goodness-of-fit measure
#' comparing the expected and observed numbers of exceedances of the
#' T-year levels for a set of return periods.
#'
#' @param xdat Numeric vector of data.
#' @param vecT Vector of return periods.
#' @param para Parameter vector.
#' @param model Model name string.
#' @return Sum of the relative absolute differences.
#' @keywords internal
gof.ene_all = function(xdat, vecT=c(5,10,20,40,80),
                       para=NULL, model=NULL){

  ns=length(xdat)
  nT = length(vecT)
  chi=rep(NA,nT)

  for(i in 1:nT){
    qt = quagev.NS(f=1-(1/vecT[i]), para,
                   nsample=ns, model)
    ene = ns/vecT
    sne = sum(xdat >= qt)
    chi[i] = abs(ene[i]-sne) /ene[i]
  }
  sum(chi)
}

#-------------------------------------------------------------
#' Modified L-moments computation with failure handling
#'
#' @description A slightly modified version of \code{lmomco::lmoms} (by
#' J. Park) which returns an \code{ifail} flag instead of stopping when
#' the L-moments cannot be computed.
#'
#' @param x Numeric vector.
#' @param nmom Number of L-moments.
#' @param mtrim If TRUE, uses left-trimmed TL-moments (leftrim=5).
#' @param no.stop If TRUE, returns \code{ifail=1} instead of stopping on error.
#' @return A list of L-moments with an additional \code{ifail} element.
#' @keywords internal
lmoms.md.park = function (x, nmom = 5, mtrim=FALSE,
                          no.stop = FALSE){
    z=list()
    n <- length(x)
    if (nmom > n) {
      if (no.stop) { z$ifail=1; return(z) }
      stop("More L-moments requested by parameter 'nmom' than data points available in 'x'")
    }
    if (length(unique(x)) == 1) {
      if (no.stop) { z$ifail=1; return(z) }
      stop("all values are equal--Lmoments can not be computed")
    }

    if(mtrim==FALSE){
      z <- TLmoms(x, nmom = nmom)
    }else{
      z <- TLmoms(x, nmom = nmom, leftrim=5)
    }
    z$source <- "lmoms"
    z$ifail=0
    return(z)
}

#---------------------------------------------------------------
#' Transform non-stationary GEV data to standard Gumbel
#'
#' @description Transforms observations from a non-stationary GEV11 model
#' with parameters \code{para} = (mu0, mu1, sigma0, sigma1, xi) to the
#' standard Gumbel scale.
#'
#' @param xdat Numeric vector of data.
#' @param para GEV11 parameter vector (mu0, mu1, sigma0, sigma1, xi).
#' @param na.ok If TRUE, invalid (out-of-support) values are returned as NA;
#'   if FALSE, they are clipped.
#' @return Numeric vector on the standard Gumbel scale.
#' @keywords internal
trans.gum01= function(xdat, para=NULL, na.ok=FALSE){

  ns= length(xdat)
  ngum01= rep(0,ns)

  ngum.dat= (xdat-(para[1]+ para[2]*seq(1,ns)))/exp(para[3]
                          + para[4]*seq(1,ns))
  work= 1 - para[5]*ngum.dat

  if(na.ok == FALSE){
    work[which(work <= 0)] = 1e-100
    ngum01= log(work)/(-para[5])
  }else{
    ngum01[which(work <= 0)] = NA
    id= which(work > 0)
    ngum01[id] = log(work[id])/(-para[5])
  }
  return(ngum01)
}
