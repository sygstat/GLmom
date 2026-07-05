# Penalty (preference) functions for the GEV shape parameter
# Shared by the stationary (glme.gev) and non-stationary (glme.gev11) methods.

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
#' @seealso \code{\link{pk.beta}} for the adaptive beta penalty,
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


#' Normal preference function for shape parameter
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
#' @seealso \code{\link{pk.beta}} for the beta penalty,
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


#' Beta preference function for the GEV shape parameter
#'
#' @description
#' Computes the data-adaptive Beta preference (penalty) function for the GEV
#' shape parameter, used by both the stationary (\code{\link{glme.gev}}) and
#' the non-stationary (\code{\link{glme.gev11}}) GLME methods (v2.0.0 unified
#' implementation).
#'
#' When \code{q} is \code{NULL} (adaptive mode), the support is
#' \eqn{(\max(-1, \hat\xi - c_0),\; \min(0.3, \hat\xi + c_0))} where
#' \eqn{\hat\xi} is \code{lme.center[3]}, and
#' \eqn{q = p + \min(|\hat\xi| c_1, c_2)} for \eqn{\hat\xi \le -0.05}
#' (otherwise \eqn{q = p}). When \code{q} is given (fixed mode), the support
#' is \code{(min.xi, max.xi)}; this mode is used for the fixed literature
#' penalties ("ms", "park", "cannon").
#'
#' \code{pk.beta.stnary} is kept as a backward-compatible alias. Note that
#' since v2.0.0 both names use this unified implementation (the v1.x default
#' hyperparameters \code{c0=0.3, c2=5} changed to \code{c0=0.35, c2=2}, and
#' out-of-support values now return \code{1e-100} instead of \code{1}).
#'
#' @param para A vector of GEV parameters; only \code{para[3]} (shape) is used.
#' @param lme.center L-moment estimates used as the center (adaptive mode).
#' @param p First shape parameter of the beta distribution.
#' @param q Second shape parameter (optional; if provided, fixed mode is used).
#' @param min.xi Lower support limit in fixed mode (default -1).
#' @param max.xi Upper support limit in fixed mode (default 0.5).
#' @param c0 Half-width of the adaptive support (default 0.35).
#' @param c1 Scaling hyperparameter for q (default 10).
#' @param c2 Upper limit hyperparameter for q (default 2).
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
#'   \code{\link{glme.gev}}, \code{\link{glme.gev11}}.
#'
#' @examples
#' # Adaptive beta preference for xi = -0.2 centered at LME xi = -0.15
#' pk.beta(para = c(100, 20, -0.2),
#'         lme.center = c(100, 20, -0.15), p = 6, c1 = 3, c2 = 1)
#'
#' @author Yonggwan Shin, Seokkap Ko, Jihong Park, Yire Shin, Jeong-Soo Park
#' @export
pk.beta = function(para=NULL, lme.center=NULL, p=NULL,
                   q=NULL, min.xi=-1.0, max.xi=0.5,
                   c0=0.35, c1=10, c2=2){

  pk=list()

  if(is.null(q)){

    aa= max(-1, lme.center[3]-c0)
    bb= min(0.3, lme.center[3]+c0)
    al=min(aa,bb)
    bl=max(aa,bb)

    if(lme.center[3] <= -0.05) {
      cut= min(lme.center[3], -0.05)
      qlim= min(abs(cut)*c1, c2)
    }else{ qlim =0.0 }

    q=p+qlim

  }else if(!is.null(q)){
    lme.center= c(0,0,0)
    al= min.xi; bl=max.xi
  }

  Be  <- integrate(Befun, lower=al, upper=bl,
                   al=al, bl=bl, p=p, q=q)[1]$value

  pk$pk.one=1e-100
  if(lme.center[3] < 0.5){
    if( (para[3] > al) & (para[3] < bl) ) {
      pk$pk.one <- ((-al+para[3])^(p-1))*((bl-para[3])^(q-1))/ Be
    }else{
      pk$pk.one = 1e-100
    }
  }
  if(is.na(pk$pk.one)) pk$pk.one=1e-100
  pk$p=p
  pk$q=q
  return(pk)
}

#' @rdname pk.beta
#' @export
pk.beta.stnary = pk.beta


#' Validate and resolve penalty settings
#'
#' @description
#' Normalizes the penalty name (e.g., "Beta" to "beta") and, when
#' \code{pen.choice} is given, resolves the preset hyperparameters:
#' for \code{pen="beta"}, choices 1-6 give (p, c1, c2) =
#' (6,3,1), (6,5,2), (6,7,3), (2,3,0.5), (2,5,1), (2,7,1.5);
#' for \code{pen="norm"}, choices 1-4 give (mu, std) =
#' (-0.5,0.25), (-0.5,0.15), (-0.6,0.25), (-0.6,0.15).
#' These presets follow Table 1 and Eq. (15) of Shin et al. (2025).
#'
#' @param pen Penalty name (character).
#' @param pen.choice Preset number or NULL.
#' @param p,c1,c2 Beta penalty hyperparameters (used when pen.choice is NULL).
#' @param mu,std Normal penalty hyperparameters (used when pen.choice is NULL).
#' @return A list with elements pen, p, c1, c2, mu, std.
#' @keywords internal
check.penalty = function(pen=NULL,pen.choice=NULL,p=NULL,
                         c1=NULL,c2=NULL,mu=NULL,std=NULL){

  if(is.null(pen) | !is.character(pen)) {
    stop("pen should be given as a character")}
  if(pen=="Beta" ) pen="beta"  ;  if(pen=="CD") pen="cd"
  if(pen=="MS") pen="ms"       ;  if(pen=="Park") pen="park"
  if(pen=="normal") pen="norm" ;  if(pen=="Cannon") pen="cannon"

  if(pen=='beta' & is.null(pen.choice)){
    if(is.null(p) | is.null(c1) | is.null(c2)){
      stop("if pen.choice is null, p, c1, c2 should be specified")}
  }
  if(pen=='norm' & is.null(pen.choice)){
    if(is.null(mu) | is.null(std)){
      stop("if pen.choice is null, mu, std should be specified")}
  }

  if(pen=='beta' & !is.null(pen.choice)){
    if(pen.choice %% 1 != 0 | pen.choice < 1 | pen.choice > 6 ){
      stop("pen.choice for beta should be an integer: 1~6")
    }
    pc1c2= matrix(c(6,6,6, 2,2,2, 3,5,7, 3,5,7,
                    1,2,3, 0.5,1,1.5), 6,3, byrow=FALSE)
    p= pc1c2[pen.choice,1]
    c1=pc1c2[pen.choice,2]
    c2=pc1c2[pen.choice,3]
  }
  if(pen=='norm' & !is.null(pen.choice)){
    if(pen.choice %% 1 != 0 | pen.choice < 1 | pen.choice > 4 ){
      stop("pen.choice for norm should be an integer: 1~4")
    }
    mustd= matrix(c(-.5,.25, -.5,.15,
                    -.6,.25, -.6,.15), 4,2, byrow=TRUE)
    mu= mustd[pen.choice,1]
    std=mustd[pen.choice,2]
  }
  result=list(pen=pen, p=p, c1=c1, c2=c2, mu=mu, std=std)
  return(result)
}


#' Penalty dispatcher (negative log prior)
#'
#' @description
#' Returns the negative log of the preference (prior) function for the
#' selected penalty type. Dispatches to \code{\link{pk.norm.stnary}}
#' ("norm"), the adaptive \code{\link{pk.beta}} ("beta"), the fixed
#' literature beta penalties ("ms": Beta(6,9); "park": Beta(2.5,2.5);
#' "cannon": Beta(2,3.3)), the Coles-Dixon exponential penalty ("cd"),
#' or zero ("no").
#'
#' @param par Parameter vector; only \code{par[3]} (shape) is penalized.
#' @param mu,std Hyperparameters for the normal penalty.
#' @param lme Center (L-moment pre-estimate) for the adaptive beta penalty.
#' @param pen Penalty type.
#' @param p,c0,c1,c2,q Hyperparameters for the beta penalty.
#' @return Negative log prior value (scalar).
#' @keywords internal
penalty.fun = function(par, mu=NULL, std=NULL, lme=NULL, pen='beta',
                       p=NULL, c0=0.35, c1=NULL, c2=NULL, q=NULL){

  if(pen=='norm' | pen=="normal"){

    pk = -log( pk.norm.stnary(para=par, mu=mu, std=std) )

  }else if(pen=='beta' | pen=="Beta"){

    pk = -log( pk.beta(para=par, lme.center=lme, p=p, q=q,
                       c0=c0, c1=c1, c2=c2)$pk.one )

  }else if(pen=='ms' | pen=="MS"){

    pk = -log( pk.beta(para=par, p=6, q=9)$pk.one )

  }else if(pen=="park" | pen=="Park"){

    pk = -log( pk.beta(para=par, p=2.5, q=2.5)$pk.one )

  }else if(pen=="cannon" | pen=="Cannon"){

    pk = -log( pk.beta(para=par, p=2, q=3.3)$pk.one )

  }else if(pen=="cd" | pen=="CD"){

    if (par[3] >= 0) { pk = 0
    }else if (par[3] > -1 & par[3] < 0) {
      pk = (1/(1 + par[3])) - 1
    }else{ pk = 10^6 }

  }else if(pen=='no'){
    pk = 0
  }else{
    stop("unknown penalty type: ", pen)
  }
  pk
}
