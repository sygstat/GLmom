# ---------------------------------------------------------------
# Generalized L-moment estimation (GLME) for the stationary GEV
# L-me under a preference (prior) function for stationary GEV
# ---------------------------------------------------------------

#' Initialize parameters for GEV MLE estimation
#'
#' @description
#' This function initializes parameters for GEV maximum likelihood estimation.
#'
#' @param data A numeric vector of data to be fitted.
#' @param ntry Number of initial parameter sets to generate.
#'
#' @details
#' The function generates `ntry` sets of initial parameters for the GEV distribution.
#' It uses L-moment estimates as a starting point and then generates additional
#' sets of parameters using random perturbations.
#'
#' @return A matrix with `ntry` rows and 3 columns, where each row represents
#' a set of initial parameters (location, scale, shape) for the GEV distribution.
#'
#' @author Yonggwan Shin, Seokkap Ko, Jihong Park, Yire Shin, Jeong-Soo Park
#' @keywords internal
init.gevmax <-function(data=NULL, ntry=NULL){

  if(ntry < 3) ntry=3
  init <-matrix(NA, ntry, 3)

  lmom_init = lmoms(data, nmom=3)
  init[1,] <- pargev(lmom_init)$para

  sd1= max(abs(init[1,1])*0.1, 5)
  sd2= max(init[1,2]*2, 2)

  maxm2=ntry-1
  init[2:ntry,1] <- init[1,1]+rnorm(n=maxm2, mean=0, sd = sd1)
  init[2:ntry,2] <- init[1,2]+rnorm(n=maxm2, mean=4, sd = sd2)
  init[2:ntry,3] <- runif(n=maxm2, min= -0.45, max=0.4)
  init[2:ntry,2] = pmax(0.1, init[2:ntry,2])

  return(init)
}


#' Initialize starting values for GLME optimization
#'
#' @description
#' Generates multiple starting parameter sets for multi-start optimization
#' in GLME estimation, for both the stationary GEV (\code{model="gev00"})
#' and the non-stationary GEV11 (\code{model="gev11"}). Uses L-moment
#' estimates as a base and adds random perturbations.
#'
#' For \code{model="gev11"} the second column is returned on the log scale
#' (matching sigma0 of the GEV11 parameterization), and when \code{pretheta}
#' is supplied the first candidate is a perturbation of the pre-estimate.
#'
#' @param data A numeric vector of data to be fitted.
#' @param ntry Number of initial parameter sets to generate (minimum 2).
#' @param model Either "gev00" (stationary, default) or "gev11".
#' @param pretheta Optional pre-estimate (mu0, mu1, sigma0, sigma1, xi)
#'   used for the first candidate when \code{model="gev11"}.
#' @param xdat Deprecated alias of \code{data} (v1.x compatibility).
#'
#' @return A matrix with \code{ntry} rows and 3 columns,
#'   where each row is a candidate starting point for optimization.
#'
#' @seealso \code{\link{glme.gev}}, \code{\link{glme.gev11}} which use this
#'   function internally.
#'
#' @examples
#' data(haenam)
#' inits <- init.glme(haenam$X1, ntry = 5)
#' print(inits)
#'
#' @author Yonggwan Shin, Seokkap Ko, Jihong Park, Yire Shin, Jeong-Soo Park
#' @export
init.glme = function(data = NULL, ntry = NULL, model="gev00",
                     pretheta=NULL, xdat=NULL)
{
  if(is.null(data) & !is.null(xdat)) data <- xdat  # v1.x compatibility

  if (ntry < 2) ntry <- 2
  init <- matrix(NA, ntry, 3)

  init[2, 1:3] <- pargev(lmoms(data, nmom = 3))$para[1:3]
  init[1, 1:3] = c(mean(data)*0.75, sd(data)*0.5,
                   max(init[2,3]-0.03,-0.99))

  if(ntry >= 3){
    maxm2 <- ntry - 2
    init[3:ntry, 1] <- init[2, 1] + rnorm(n = maxm2, mean = 0,
                                       sd = 5)
    init[3:ntry, 2] <- pmax(0.1, init[2, 2] + rnorm(n = maxm2,
                                      mean = 3, sd = 3))
    init[3:ntry, 3] <- runif(n = maxm2, min = -0.49, max = 0.4)
  }

  if(model=="gev11"){
    init[,2] <- log(init[,2])

    if(!is.null(pretheta)){
      init[1, 1:3] = c( pretheta[1]*0.98, max(pretheta[3]-0.1, 0.1),
                        max(pretheta[5]-0.03, -0.99) )
    }
  }
  return(init)
}


#' Bootstrap covariance of sample L-moments
#'
#' @description Computes the covariance matrix of the first three sample
#' L-moments by nonparametric bootstrap. Used as a fallback when the
#' direct covariance estimate is singular.
#'
#' @param xdat Numeric data vector.
#' @param BB Number of bootstrap replicates.
#' @return A 3x3 covariance matrix.
#' @keywords internal
boot.cov = function(xdat,BB){
  sam.lmom= matrix(NA,BB,3)
  ns= length(xdat)

  for(ib in 1:BB){
    idx= sample.int(ns,ns,replace=TRUE)
    sam.lmom[ib,1:3]= lmoms(xdat[idx],nmom=3)$lambdas
  }
  cov(sam.lmom)
}


#' Calculate the likelihood for Generalized L-moments estimation of GEV distribution
#'
#' @description
#' This function calculates the likelihood (or more precisely, a penalized negative log-likelihood)
#' for the Generalized L-moments estimation of the Generalized Extreme Value (GEV) distribution.
#'
#' @param par A vector of GEV parameters (location, scale, shape).
#' @param xdat A numeric vector of data.
#' @param slmgev Sample L-moments of the data.
#' @param covinv Inverse of the covariance matrix of the sample L-moments.
#' @param lcovdet Log determinant of the covariance matrix.
#' @param mu Mean for the normal penalization (used when pen='norm').
#' @param std Standard deviation for the normal penalization (used when pen='norm').
#' @param lme L-moment estimates of the parameters.
#' @param pen Penalization method: 'norm', 'beta', 'ms', 'park', 'cannon', 'cd', or 'no'.
#' @param p Shape parameter for beta penalty.
#' @param c1 Scaling parameter for beta penalty.
#' @param c2 Upper limit parameter for beta penalty.
#' @param c0 Half-width of the adaptive beta penalty support (default 0.35).
#' @param q Optional fixed second shape parameter for the beta penalty.
#'
#' @details
#' The function performs the following steps:
#' 1. Checks if the parameters are within valid ranges.
#' 2. Calculates the expected L-moments based on the current parameters.
#' 3. Computes the difference between expected and sample L-moments.
#' 4. Calculates the generalized L-moments distance.
#' 5. Applies a penalization term based on the specified method.
#' 6. Returns the sum of the L-moments distance and the penalization term.
#'
#' @return A numeric value representing the penalized negative log-likelihood.
#' A lower value indicates a better fit.
#'
#' @references
#' Shin, Y., Shin, Y., Park, J. & Park, J.-S. (2025). Generalized method of
#' L-moment estimation for stationary and nonstationary extreme value models.
#' arXiv preprint arXiv:2512.20385. \doi{10.48550/arXiv.2512.20385}
#'
#' @seealso \code{\link{glme.gev}} which calls this function for optimization.
#'
#' @examples
#' data(haenam)
#' x <- haenam$X1
#' slm <- lmomco::lmoms(x, nmom = 3)
#' cov_mat <- lmomco::lmoms.cov(x, nmom = 3)
#' lme_par <- lmomco::pargev(slm)$para
#' glme.like(par = lme_par, xdat = x, slmgev = slm,
#'           covinv = solve(cov_mat), lcovdet = log(det(cov_mat)),
#'           mu = -0.5, std = 0.2, lme = lme_par, pen = "beta",
#'           p = 6, c1 = 3, c2 = 1)
#'
#' @author Yonggwan Shin, Seokkap Ko, Jihong Park, Yire Shin, Jeong-Soo Park
#' @export
glme.like = function(par, xdat=xdat, slmgev=slmgev, covinv=covinv,
                     lcovdet=lcovdet, mu=mu, std=std, lme=lme,
                     pen=pen, p=p, c1=c1, c2=c2, c0=0.35, q=NULL){

  if( par[2] <= 0) return(10^8)
  if( abs(par[3]) >= 1) return(10^8)
  if( abs(par[3]) <= 1e-5) par[3]= -1e-5

  emom= lmomgev( vec2par(par,type='gev') )$lambdas[1:3]
  if( any(is.na(emom[1:3])) ) return(10^8)

  zvec= emom - slmgev$lambdas[1:3]

  gld= t(zvec) %*% covinv %*% zvec
  prob.norm =  gld/2  + (3/2)*log(2*pi) + lcovdet

  pk_beta = penalty.fun(par=par, mu=mu, std=std, lme=lme,
                    pen=pen, p=p, c0=c0, c1=c1, c2=c2, q=q)
  as.numeric(prob.norm + pk_beta)
}


#' Generalized L-moments estimation for generalized extreme value distribution
#'
#' @description
#' This function estimates the parameters of the Generalized Extreme Value
#' distribution by the generalized method of L-moment estimation (GLME),
#' which minimizes the generalized L-moment distance plus a penalty
#' (preference) function on the shape parameter.
#'
#' @param xdat A numeric vector of data to be fitted.
#' @param ntry Number of attempts for parameter estimation (default 5). Higher
#'   values increase the chance of finding the global optimum by trying
#'   different initial conditions.
#' @param pen Type of penalty function: Choose among "norm", "beta" (default),
#'   "ms", "park", "cannon", "cd", and "no" (without penalty function).
#' @param pen.choice Choice number of penalty function specifying hyperparameters
#'   (default 1). For "beta": 1-6 correspond to (p, c1, c2) =
#'   (6,3,1), (6,5,2), (6,7,3), (2,3,0.5), (2,5,1), (2,7,1.5).
#'   For "norm": 1-4 correspond to (mu, std) =
#'   (-0.5,0.25), (-0.5,0.15), (-0.6,0.25), (-0.6,0.15).
#'   Set \code{pen.choice=NULL} to use the hyperparameters given by
#'   \code{p, c1, c2} (beta) or \code{mu, std} (norm) directly.
#' @param mu Mean hyperparameter for "norm" penalty function (default -0.5).
#' @param std Standard deviation hyperparameter for "norm" penalty function (default 0.2).
#' @param p Shape hyperparameter for "beta" penalty function (default 6).
#' @param c1 Scaling hyperparameter for "beta" penalty function (default 3).
#' @param c2 Upper limit hyperparameter for "beta" penalty function (default 1).
#' @param c0 Half-width of the adaptive beta penalty support (default 0.35).
#' @param q Optional fixed second shape parameter for the beta penalty. If
#'   given, the beta penalty support is fixed instead of data-adaptive.
#' @param show If TRUE, prints the objective value and parameters of each try.
#' @param method Optimization method passed to \code{\link[stats]{optim}}
#'   (default "BFGS").
#' @param maxit Maximum number of iterations for \code{optim} (default 70).
#' @param abstol Absolute and relative convergence tolerance for \code{optim}
#'   (default 1e-5).
#'
#' @details
#' The equations for the L-moments for LME of the GEVD are
#' \deqn{ \underline{\bf \lambda} - \underline{\bf l} = \underline{\bf 0},}
#' where \eqn{ \underline{\bf \lambda} =(\lambda_1,\; \lambda_2,\; \lambda_3)^t } and \eqn{\underline{\bf l} =(l_1,\; l_2,\; l_3)^t}.
#' Next, we define the generalized L-moments distance (GLD) as;
#' \deqn{(\underline{\bf \lambda} -\underline{\bf l})^t V^{-1} (\underline{\bf \lambda} -\underline{\bf l}),}
#' where \eqn{V} is the variance-covariance matrix of the sample L-moments up to the third order.
#'
#' @return The glme.gev function returns a list containing the following elements:
#' \itemize{
#'  \item para.glme - The estimated parameters of the Generalized Extreme Value distribution.
#'  \item para.lme - The L-moment estimates of the parameters.
#'  \item nllh.glme - The penalized negative log-likelihood of the GLME solution.
#'  \item convergence - 0 if the optimization converged, 5 otherwise.
#'  \item pen, pen_pen.choice - The penalization method (and choice) used.
#'  \item p_q - (for beta penalty) The p and q values used.
#'  \item c1_c2, c0_c1_c2 - (for beta penalty) The hyperparameters used.
#'  \item mu_std - (for norm penalty) The mu and std values used.
#'  \item covinv.lmom - The inverse of the covariance matrix of the L-moments.
#'  \item lcovdet - The log determinant of the covariance matrix.
#'  \item data - The input data (used by the plot method).
#' }
#' The object has class \code{"glme"}; see \code{\link{GLmom-methods}}.
#'
#' @references
#' Shin, Y., Shin, Y., Park, J. & Park, J.-S. (2025). Generalized method of
#' L-moment estimation for stationary and nonstationary extreme value models.
#' arXiv preprint arXiv:2512.20385. \doi{10.48550/arXiv.2512.20385}
#'
#' @seealso \code{\link{glme.gev11}} for non-stationary GEV estimation,
#'   \code{\link{ma.gev}} for model averaging estimation,
#'   \code{\link{glme.like}} for the objective function,
#'   \code{\link{quagev.NS}} for quantile computation.
#'
#' @author Yonggwan Shin, Seokkap Ko, Jihong Park, Yire Shin, Jeong-Soo Park
#'
#' @examples
#' # Load example heavy-rainfall data
#' data(haenam)
#' x <- haenam$X1
#'
#' # Estimate GEV parameters using beta penalty (default)
#' result <- glme.gev(x, ntry = 5)
#' print(result$para.glme)
#'
#' # Using Martins-Stedinger penalty
#' result_ms <- glme.gev(x, ntry = 5, pen = "ms")
#' print(result_ms$para.glme)
#'
#' @export
glme.gev= function(xdat, ntry=5, pen='beta', pen.choice=1,
                   mu= -0.5, std= 0.2, p=6, c1=3, c2=1,
                   c0=0.35, q=NULL, show=FALSE,
                   method="BFGS", maxit=70, abstol=1e-5){
  # updated on 17Apr26 (v2.0.0, revised code of Shin et al.)

  z=list()

  ch= check.penalty(pen=pen, pen.choice=pen.choice, p=p,
                   c1=c1,c2=c2,mu=mu,std=std)
  pen= ch$pen
  if(!is.null(pen.choice)){
    p= ch$p ; c1 =ch$c1;  c2= ch$c2
    mu= ch$mu ; std= ch$std
  }

  sinit <- init.glme(xdat, ntry=ntry)
  slmgev = lmoms(xdat, nmom=3)

  z$para.lme= lme= pargev(slmgev)$para
  conv= nllh= rep(1e6, ntry)
  para= matrix(NA,ntry,3)

  cov= lmoms.cov(xdat, nmom=3)
  detc = det(cov)

  if(detc <= 0){    # Bootstrap to calculate cov
    cov= boot.cov(xdat,BB=100)
    detc= det(cov)
  }
  lcovdet= log(detc)
  covinv = solve(cov)

  z$covinv.lmom = covinv   # kept for v1.x compatibility
  z$lcovdet = lcovdet

  # estimating paras using optim
  for(i in 1:ntry){

    kt <- try(
      optim(par=as.vector(sinit[i,1:3]), fn=glme.like,
            method=method,
            control=list(abstol=abstol, reltol=abstol,
                         maxit=maxit),
            xdat=xdat, slmgev=slmgev, covinv=covinv,
            lcovdet=lcovdet, mu=mu, std=std, lme=lme,
            pen=pen, p=p, c1=c1, c2=c2, c0=c0, q=q),
      silent=TRUE)

    if(inherits(kt, "try-error")){
      nllh[i]= 1e6
      conv[i]= 5
      next
    }

    if( kt$convergence != 0){
      nllh[i]=10^6
      conv[i]=5
    }else{
      nllh[i]= kt$value
      conv[i]= 0
      para[i,]= kt$par
    }

    if(show==TRUE){
      cat("--- itry, nllh=",i, nllh[i],"\n")
      cat("--- itry, para=",i, round(kt$par,4),"\n","\n")
    }
  } #for

  if(all(nllh >= 10^5)) {
    message("-- No solution found in optim at glme.gev")
    z$para.glme = z$para.lme
    z$convergence = 5
    z$pen = pen
    z$data = xdat
    class(z) = "glme"
    return(z)
  }
  istar= which.min(nllh)
  z$para.glme = para[istar,]
  z$nllh.glme = nllh[istar]
  z$convergence = conv[istar]

  names(z$para.lme) <- names(z$para.glme) <- c("mu","sig","xi")
  z$pen = pen                          # kept for v1.x compatibility
  z$pen_pen.choice = c(pen, pen.choice)

  if(pen=="beta"){
    ww = pk.beta(para= z$para.glme, lme.center=lme, p=p,
                 c0=c0, c1=c1, c2=c2, q=q)
    z$p_q= c(ww$p, ww$q)
    z$c1_c2=c(c1,c2)                   # kept for v1.x compatibility
    z$c0_c1_c2=c(c0,c1,c2)

  }else if(pen=="norm"){
    z$mu_std = c(mu, std)
  }
  z$data = xdat
  class(z) = "glme"
  return(z)
}


#' GEV parameter estimation with fixed shape parameter
#'
#' @description
#' Estimates GEV location and scale parameters from L-moments while keeping
#' the shape parameter fixed at a user-specified value. Modified from
#' \code{lmomco::pargev()}.
#'
#' @param lmom L-moments object.
#' @param kfix Fixed shape parameter value.
#' @param checklmom Whether to check L-moment validity.
#' @param ... Additional arguments.
#' @return A list with components:
#' \describe{
#'   \item{type}{Character "gev"}
#'   \item{para}{Numeric vector of GEV parameters (xi=location, alpha=scale, kappa=shape)}
#'   \item{source}{Character "pargev"}
#' }
#'
#' @references
#' Hosking, J. R. M. (1990). L-moments: Analysis and estimation of
#' distributions using linear combinations of order statistics.
#' Journal of the Royal Statistical Society, Series B, 52(1), 105-124.
#' \doi{10.1111/j.2517-6161.1990.tb01775.x}
#'
#' @seealso \code{\link{glme.gev}} for GLME estimation,
#'   \code{\link[lmomco]{pargev}} for the original L-moment GEV fitting.
#'
#' @examples
#' data(haenam)
#' lmom <- lmomco::lmoms(haenam$X1, nmom = 3)
#' pargev.kfix(lmom, kfix = -0.1)
#'
#' @author Yonggwan Shin, Seokkap Ko, Jihong Park, Yire Shin, Jeong-Soo Park
#' @export
pargev.kfix= function (lmom, kfix= 0.1, checklmom = TRUE, ...)
{

  # modified from 'pargev' function in lmomco package

  para <- rep(NA, 3)
  names(para) <- c("xi", "alpha", "kappa")
  SMALL <- 1e-05
  EPS <- 1e-06
  MAXIT <- 20
  EU <- 0.57721566
  DL2 <- 0.69314718
  DL3 <- 1.0986123
  A0 <- 0.2837753
  A1 <- -1.21096399
  A2 <- -2.50728214
  A3 <- -1.13455566
  A4 <- -0.07138022
  B1 <- 2.06189696
  B2 <- 1.31912239
  B3 <- 0.25077104
  C1 <- 1.59921491
  C2 <- -0.48832213
  C3 <- 0.01573152
  D1 <- -0.64363929
  D2 <- 0.08985247
  if (length(lmom$L1) == 0) {
    lmom <- lmorph(lmom)
  }
  if (checklmom & !are.lmom.valid(lmom)) {
    warning("L-moments are invalid")
    return()
  }
  T3 <- lmom$TAU3
  if (T3 > 0) {
    Z <- 1 - T3
    G <- (-1 + Z * (C1 + Z * (C2 + Z * C3)))/(1 + Z * (D1 +
                                                         Z * D2))
    if (abs(G) < SMALL) {
      para[3] <- 0
      para[2] <- lmom$L2/DL2
      para[1] <- lmom$L1 - EU * para[2]
      return(list(type = "gev", para = para))
    }
  }
  else {
    G <- (A0 + T3 * (A1 + T3 * (A2 + T3 * (A3 + T3 * A4))))/(1 +
                                                               T3 * (B1 + T3 * (B2 + T3 * B3)))
    if (T3 >= -0.8) {
    }
    else {
      if (T3 <= -0.97)
        G <- 1 - log(1 + T3)/DL2
      T0 <- (T3 + 3) * 0.5
      CONVERGE <- FALSE
      for (it in seq(1, MAXIT)) {
        X2 <- 2^-G
        X3 <- 3^-G
        XX2 <- 1 - X2
        XX3 <- 1 - X3
        T <- XX3/XX2
        DERIV <- (XX2 * X3 * DL3 - XX3 * X2 * DL2)/(XX2 *
                                                      XX2)
        GOLD <- G
        G <- G - (T - T0)/DERIV
        if (abs(G - GOLD) <= EPS * G)
          CONVERGE <- TRUE
      }
      if (CONVERGE == FALSE) {
        warning("Noconvergence---results might be unreliable")
      }
    }
  }

  para[3]=kfix
  G = kfix

  GAM <- exp(lgamma(1 + G))
  para[2] <- lmom$L2 * G/(GAM * (1 - 2^(-G)))
  para[1] <- lmom$L1 - para[2] * (1 - GAM)/G
  return(list(type = "gev", para = para, source = "pargev"))
}
