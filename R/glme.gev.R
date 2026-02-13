# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


#----------------------------------------------------------------
# L-me under a preference (prior) function for stationary GEV ---
# ---------------------------------------------------------------

#' Beta function integrand for penalty calculation
#'
#' @param w Integration variable.
#' @param al Lower bound.
#' @param bl Upper bound.
#' @param p Shape parameter p.
#' @param q Shape parameter q.
#' @return Integrand value.
#' @keywords internal
Befun <- function(w, al=al, bl=bl, p=p, q=q) {
  ((-al+w)^(p-1)) * ((bl-w)^(q-1))
}


#' Martins-Stedinger prior function
#'
#' @description
#' Computes the Martins-Stedinger Beta(6,9) prior probability for the GEV
#' shape parameter on the interval \eqn{[-0.5, 0.5]}.
#'
#' @param para A vector of GEV parameters (location, scale, shape).
#' @param p Shape parameter for beta distribution (default 6).
#' @param q Shape parameter for beta distribution (default 9).
#' @return Prior probability value (scalar).
#'
#' @references
#' Martins, E. S. & Stedinger, J. R. (2000). Generalized maximum-likelihood
#' generalized extreme-value quantile estimators for hydrologic data.
#' Water Resources Research, 36(3), 737-744.
#' \doi{10.1029/1999WR900330}
#'
#' @seealso \code{\link{pk.beta.stnary}} for the adaptive beta penalty,
#'   \code{\link{glme.gev}} which uses these penalty functions.
#'
#' @examples
#' # Evaluate MS prior at xi = -0.2
#' MS_pk(para = c(100, 20, -0.2))
#'
#' @author Yonggwan Shin, Seokkap Ko, Jihong Park, Yire Shin, Jeong-Soo Park
#' @export
MS_pk = function(para=para, p=6, q=9){

  Bef <- function(x) { ((0.5+x)^(p-1)) * ((0.5-x)^(q-1)) }
  Be  <- integrate(Bef, lower=-0.5, upper=0.5)[1]$value

  if( abs(para[3]) < 0.5 ){
    pk.ms <- ((0.5+para[3])^(p-1))*((0.5-para[3])^(q-1))/ Be
  }else if ( abs(para[3]) >= 0.5 ) {
    pk.ms = 1e-20
  }

  pk.ms
}


#' Normal preference function for shape parameter (stationary GEV)
#'
#' @description
#' Computes a normal distribution-based preference (penalty) function value
#' for the GEV shape parameter. The alias \code{new_pf_norm} is provided
#' for backward compatibility.
#'
#' @param para A vector of GEV parameters.
#' @param mu Mean for normal distribution.
#' @param std Standard deviation for normal distribution.
#' @return Preference function value (scalar).
#'
#' @references
#' Shin, Y., Shin, Y., Park, J. & Park, J.-S. (2025). Generalized method of
#' L-moment estimation for stationary and nonstationary extreme value models.
#' arXiv preprint arXiv:2512.20385. \doi{10.48550/arXiv.2512.20385}
#'
#' @seealso \code{\link{pk.beta.stnary}} for the beta penalty,
#'   \code{\link{MS_pk}} for the Martins-Stedinger penalty,
#'   \code{\link{glme.gev}} which uses these penalty functions.
#'
#' @examples
#' # Normal preference with mean=-0.5, sd=0.2 at xi=-0.2
#' pk.norm.stnary(para = c(100, 20, -0.2), mu = -0.5, std = 0.2)
#'
#' @author Yonggwan Shin, Seokkap Ko, Jihong Park, Yire Shin, Jeong-Soo Park
#' @export
pk.norm.stnary = function(para=NULL, mu=NULL, std=NULL){
  1 + dnorm(para[3], mean= mu, sd=std)
}

#' @rdname pk.norm.stnary
#' @export
new_pf_norm = pk.norm.stnary


#' Beta preference function for stationary GEV
#'
#' @description
#' Computes a Beta distribution-based adaptive preference (penalty) function
#' for the GEV shape parameter. The hyperparameters are adapted based on the
#' L-moment estimate of the shape parameter.
#'
#' @param para A vector of GEV parameters.
#' @param lme.center L-moment estimates as center.
#' @param p Shape parameter (default 6).
#' @param q Shape parameter q (optional, if provided uses fixed limits).
#' @param c0 Limit parameter (default 0.3).
#' @param c1 Scaling parameter (default 10).
#' @param c2 Upper limit parameter (default 5).
#' @return A list containing:
#' \describe{
#'   \item{pk.one}{Preference function value (scalar)}
#'   \item{p}{Shape parameter p used}
#'   \item{q}{Shape parameter q used}
#' }
#'
#' @references
#' Shin, Y., Shin, Y., Park, J. & Park, J.-S. (2025). Generalized method of
#' L-moment estimation for stationary and nonstationary extreme value models.
#' arXiv preprint arXiv:2512.20385. \doi{10.48550/arXiv.2512.20385}
#'
#' @seealso \code{\link{pk.norm.stnary}} for the normal penalty,
#'   \code{\link{MS_pk}} for the Martins-Stedinger penalty,
#'   \code{\link{glme.gev}} which uses these penalty functions.
#'
#' @examples
#' # Beta preference for xi = -0.2 centered at LME xi = -0.15
#' pk.beta.stnary(para = c(100, 20, -0.2),
#'                lme.center = c(100, 20, -0.15), p = 6)
#'
#' @author Yonggwan Shin, Seokkap Ko, Jihong Park, Yire Shin, Jeong-Soo Park
#' @export
pk.beta.stnary = function(para=NULL, lme.center=NULL, p=NULL,
                          q=NULL, c0=0.3, c1=10, c2=5){

  pk=list()
  pk.one = 1e-10

  if(is.null(q)){

    ulim= c0                      #min(0.3, abs(lme[3])*3 )
    aa= max(-1, lme.center[3]-ulim)
    bb= min(0.5, lme.center[3]+ulim)
    al=min(aa,bb)
    bl=max(aa,bb)

    if(lme.center[3] <= 0) {
      qlim= min(abs(lme.center[3])*c1, c2)
    }else{ qlim =0.0 }

    q=p+qlim

  }else if(!is.null(q)){
    lme.center[3]=0
    al= -0.5; bl=0.5
  }

  Be  <- integrate(Befun, lower=al, upper=bl,
                   al=al, bl=bl, p=p, q=q)[1]$value

  pk$pk.one=1
  if(lme.center[3] < 0.5){
    if( (para[3] > al) & (para[3] < bl) ) {
      pk$pk.one <- ((-al+para[3])^(p-1))*((bl-para[3])^(q-1))/ Be
    }}
  if(is.na(pk$pk.one)) pk$pk.one=1
  pk$p=p
  pk$q=q
  return(pk)
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
#' data(streamflow)
#' x <- streamflow$r1
#' slm <- lmomco::lmoms(x, nmom = 3)
#' cov_mat <- lmomco::lmoms.cov(x, nmom = 3)
#' lme_par <- lmomco::pargev(slm)$para
#' glme.like(par = lme_par, xdat = x, slmgev = slm,
#'           covinv = solve(cov_mat), lcovdet = log(det(cov_mat)),
#'           mu = -0.5, std = 0.2, lme = lme_par, pen = "beta",
#'           p = 6, c1 = 10, c2 = 5)
#'
#' @author Yonggwan Shin, Seokkap Ko, Jihong Park, Yire Shin, Jeong-Soo Park
#' @export
glme.like = function(par=par, xdat=xdat, slmgev=slmgev, covinv=covinv,
                     lcovdet=lcovdet, mu=mu, std=std, lme=lme, pen=pen,
                     p=p, c1=c1, c2=c2){

  if( par[2] <= 0) return(10^8)
  if( abs(par[3]) > 1) return(10^8)

  nsample=length(xdat)

  if( abs(par[3]) < 1e-5) par[3]= -(1e-4)

  emom= lmomgev( vec2par(par,type='gev') )

  if( is.na(emom$lambda[1]) ) return(10^8)
  if( is.na(emom$lambda[2]) ) return(10^8)
  if( is.na(emom$ratios[3]) ) return(10^8)

  zvec=rep(NA,3)
  zvec= emom$lambdas[1:3] - slmgev$lambdas[1:3]

  if( any(is.na(zvec)) ) return(10^8)

  z= t(zvec) %*% covinv %*% zvec

  prob.norm =  z/2   + (3/2)*log( (2*pi) ) + lcovdet

  if(pen=='norm' | pen=="normal"){

    pk_beta = -log( pk.norm.stnary(para=par, mu= mu, std= std) )

  }else if(pen=='beta' | pen=="Beta"){

    work = pk.beta.stnary(para= par, lme.center=lme, p=p,
                          c1=c1, c2=c2)$pk.one
    pk_beta = -log(work)

  }else if(pen=='ms' | pen=="MS"){

    pk_beta = -log( pk.beta.stnary(para=par, p=6, q=9)$pk.one)

  }else if(pen=="park" | pen=="Park"){

    pk_beta = -log( pk.beta.stnary(para=par, p=2.5, q=2.5)$pk.one)

  }else if(pen=="cannon" | pen=="Cannon"){

    pk_beta = -log( pk.beta.stnary(para=par, p=2, q=3.3)$pk.one  )

  }else if(pen=="cd" | pen=="CD"){

    if (par[3] >= 0) {pk_beta <- 0
    }else if (par[3] > -1 & par[3] < 0) {
      pk_beta <- -log( exp(-((1/(1 + par[3])) - 1)) )
    }else if (par[3] <= -1) {pk_beta = 10^6
    }

  }else if(pen=='no'){
    pk_beta =0
  }

  zz= prob.norm  + pk_beta

  return(zz)
}


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

#' Initialize random starting values for GLME optimization
#'
#' @description
#' Generates multiple random starting parameter sets for multi-start
#' optimization in GLME estimation. Uses L-moment estimates as a base
#' and adds random perturbations.
#'
#' @param xdat A numeric vector of data to be fitted.
#' @param ntry Number of initial parameter sets to generate.
#'
#' @return A matrix with \code{ntry} rows and 3 columns (location, scale, shape),
#'   where each row is a candidate starting point for optimization.
#'
#' @seealso \code{\link{glme.gev}} which uses this function internally.
#'
#' @examples
#' data(streamflow)
#' inits <- init.glme(streamflow$r1, ntry = 5)
#' print(inits)
#'
#' @author Yonggwan Shin, Seokkap Ko, Jihong Park, Yire Shin, Jeong-Soo Park
#' @export
init.glme <-function(xdat, ntry=ntry){
  init.gevmax(data=xdat, ntry=ntry)
}


#' Generalized L-moments estimation for generalized extreme value distribution
#'
#' @description
#' This function estimates the Generalized L-moments of Generalized Extreme Value distribution.
#'
#' @param xdat A numeric vector of data to be fitted.
#' @param ntry Number of attempts for parameter estimation. Higher values increase
#'   the chance of finding a more accurate estimate by trying different initial conditions.
#' @param pen Type of penalty function: Choose among "norm", "beta" (default),
#'   "ms", "park", "cannon", "cd", and "no" (without penalty function).
#' @param pen.choice Choice number of penalty function specifying hyperparameters.
#'   For "beta": 1-6 correspond to different (p, c1, c2) combinations.
#'   For "norm": 1-4 correspond to different (mu, std) combinations.
#' @param mu Mean hyperparameter for "norm" penalty function (default -0.5).
#' @param std Standard deviation hyperparameter for "norm" penalty function (default 0.2).
#' @param p Shape hyperparameter for "beta" penalty function (default 6).
#' @param c1 Scaling hyperparameter for "beta" penalty function (default 10).
#' @param c2 Upper limit hyperparameter for "beta" penalty function (default 5).
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
#'  \item covinv.lmom - The inverse of the covariance matrix of the L-moments.
#'  \item lcovdet - The log determinant of the covariance matrix.
#'  \item nllh.glme - The negative log-likelihood of the GLME solution.
#'  \item pen - The penalization method used.
#'  \item p_q - (for beta penalty) The p and q values used.
#'  \item c1_c2 - (for beta penalty) The c1 and c2 values used.
#'  \item mu_std - (for norm penalty) The mu and std values used.
#' }
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
#' # Load example streamflow data
#' data(streamflow)
#' x <- streamflow$r1
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
glme.gev= function(xdat, ntry=10, pen='beta', pen.choice=NULL,
                   mu= -0.5, std= 0.2, p=6, c1=10, c2=5){

  z=list()
  k =list()
  if(is.null(pen) | !is.character(pen)) {
    stop("pen should be given as a character")}

  if(pen=='beta' & !is.null(pen.choice)){
    if(pen.choice %% 1 != 0 | pen.choice < 1 | pen.choice > 6 ){
      stop("pen.choice for beta should be an integer: 1~6")
    }
    pc1c2= matrix(c(6,6,6,2,2,2,10,20,30,10,20,30,5,7,9,5,7,9),
                  6,3, byrow=FALSE)
    p= pc1c2[pen.choice,1]
    c1=pc1c2[pen.choice,2]
    c2=pc1c2[pen.choice,3]
  }
  if(pen=='norm' & !is.null(pen.choice)){
    if(pen.choice %% 1 != 0 | pen.choice < 1 | pen.choice > 4 ){
      stop("pen.choice for norm should be an integer: 1~4")
    }
    mustd= matrix(c(-.5,.2,-.5,.1,-.6,.2,-.6,.1),4,2, byrow=TRUE)
    mu= mustd[pen.choice,1]
    std=mustd[pen.choice,2]
  }

  # initial setting ------
  nsample=length(xdat)
  sinit=matrix(0, nrow=ntry, ncol=3)

  sinit <- init.gevmax(xdat, ntry=ntry)

  lmom_init = lmoms(xdat)
  lmom_est <- pargev(lmom_init)

  lme = lmom_est$para
  z$para.lme =  lmom_est$para

  precis=rep(NA, ntry)
  pk.ms=rep(NA, ntry)
  pk.ms.lme=rep(NA, ntry)
  isol=0
  sol=list()
  mindist=1000
  dist=rep(1000, ntry)

  covinv= matrix(NA, 3, 3)

  slmgev=lmoms(xdat)
  cov=lmoms.cov(xdat, nmom=3)

  covinv=solve(cov)
  detc = det(cov)

  #--------------------------------------------------
  if(detc <= 0){

    BB=200          # we need Bootstrap to calculate cov ---
    sam.lmom= matrix(NA,BB,3)

    for (ib in 1:BB){
      sam.lmom[ib,1:3]=lmoms(sample(xdat,size=nsample,replace=TRUE),
                             nmom=3)$lambdas
    }
    cov=cov(sam.lmom)
    covinv=solve(cov)
    detc=det(cov)
  }

  lcovdet=log(detc)
  z$covinv.lmom =covinv
  z$lcovdet =lcovdet

  #-------------------------------------------------------
  # estimating paras using nleqslv or optim
  tryCatch(
    for(i in 1:ntry){

      value=list()

      value <- try(
        optim(par=as.vector(sinit[i,1:3]), fn=glme.like,
              xdat=xdat, slmgev=slmgev, covinv=covinv,
              lcovdet=lcovdet, mu=mu, std=std, lme=lme, pen=pen,
              p=p, c1=c1, c2=c2)
      )

      if(is(value)[1]=="try-error"){
        k[[i]]$fvec <- 10^6
      }else{
        k[[i]] <- value
        k[[i]]$root = value$par
        k[[i]]$fvec = value$value
      }

      if( value$convergence != 0) {precis[i]=10^6
      }else{
        isol=isol+1
        precis[i] = k[[i]]$fvec
      }

    } #for
  ) #tryCatch

  if(isol==0) {
    message("-- No solution was found in nleqslv or optim --")
    z$para.glme = z$para.lme
    return(z)
  }

  selc_num = which.min( precis )    #precis=k[[i]]$fvec

  x  <- k[[selc_num]]

  z$para.glme = x$root
  z$nllh.glme = k[[selc_num]]$fvec
  z$pen = pen

  if(pen=="beta" | pen=="Beta"){
    ww = pk.beta.stnary(para= z$para.glme, lme.center=lme, p=p,
                        c1=c1, c2=c2)
    z$p_q= c(ww$p, ww$q); z$c1_c2=c(c1,c2)

  }else if(pen=="norm" | pen=="normal"){
    z$mu_std = c(mu, std)
  }

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
#' data(streamflow)
#' lmom <- lmomco::lmoms(streamflow$r1, nmom = 3)
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
