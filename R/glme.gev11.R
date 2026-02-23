# R code for Non-Stationary GEV11 model using L-moment
# GEV11: mu(t) = mu0 + mu1*t, sigma(t) = exp(sigma0 + sigma1*t)

#---------------------------------------------------------------------
#' Beta preference function for non-stationary GEV
#'
#' @param para A vector of GEV parameters.
#' @param lme.center L-moment estimates as center.
#' @param p Shape parameter.
#' @param c0 Limit parameter (default 0.3).
#' @param c1 Scaling parameter (default 10).
#' @param c2 Upper limit parameter (default 5).
#' @return A list containing pk.one, p, and q.
#' @author Yonggwan Shin, Seokkap Ko, Jihong Park, Yire Shin, Jeong-Soo Park
#' @keywords internal
pk.beta.ns = function(para=NULL, lme.center=NULL, p=NULL,
                      c0=0.3, c1=10, c2=5){

  pk=list()
  pk.one = 1e-10
  ulim= c0
  aa= max(-1, lme.center[3]-ulim)
  bb= min(0.3, lme.center[3]+ulim)
  al=min(aa,bb)
  bl=max(aa,bb)

  if(lme.center[3] <= 0) {
    qlim= min(0.0+abs(lme.center[3])*c1, c2)
  }else{ qlim =0.0 }

  p=p; q=p+qlim

  Be  <- integrate(Befun, lower=al, upper=bl,
                   al=al, bl=bl, p=p, q=q)[1]$value

  pk$pk.one=1e-50
  if(lme.center[3] <= 0.3){
    if( (para[3] > al) & (para[3] < bl) ) {
      pk$pk.one <- ((-al+para[3])^(p-1))*((bl-para[3])^(q-1))/ Be
    }}
  if(is.na(pk$pk.one)) pk$pk.one=1e-50

  pk$p =p; pk$q=q
  return(pk)
}

#---------------------------------------------------------------------
#' GLME objective function for GEV11 model (mu0, sigma0, xi optimization)
#'
#' @param a Parameter vector (mu0, sigma0, xi).
#' @param xdat Data vector.
#' @param newtheta Full parameter vector (mu0, mu1, sigma0, sigma1, xi).
#' @param covinv Inverse covariance matrix.
#' @param lcovdet Log determinant of covariance.
#' @param pen Penalty type.
#' @param mu Normal penalty mean.
#' @param std Normal penalty std.
#' @param p Beta penalty shape.
#' @param c1 Beta penalty scaling.
#' @param c2 Beta penalty limit.
#' @return Negative log-likelihood value.
#' @author Yonggwan Shin, Seokkap Ko, Jihong Park, Yire Shin, Jeong-Soo Park
#' @keywords internal
nllh.glme.gev11 <- function(a, xdat=xdat, newtheta=newtheta,
                            covinv=covinv, lcovdet=lcovdet,
                            pen=pen, mu=mu, std=std,
                            p=p, c1=c1, c2=c2)
{
  zvec=rep(100,3)

  mu0 <- a[1]
  mu1 <- newtheta[2]
  sig0 <- a[2]
  sig1 = newtheta[4]
  xi= a[3]

  lme.center= c(newtheta[1], newtheta[3], newtheta[5])

  ns=length(xdat)
  year=seq(1,ns)
  gum01=rep(NA, ns)
  gum.dat=rep(NA, ns)
  newg2=rep(NA,ns)

  gum.dat[1:ns]= xdat[1:ns]-(mu0 + mu1*year[1:ns])
  gum.dat[1:ns]= gum.dat[1:ns]/exp(sig0 + sig1*year[1:ns])

  gum01[1:ns]= 1-xi*gum.dat[1:ns]

  for (it in 1:ns ) {
    if( is.na(gum01[it]) ){
      newg2[it]=NA
    }else if( gum01[it] <= 0 ) {
      newg2[it]= NA
    }else if( gum01[it] > 0) {
      newg2[it]= log(gum01[it])/(-xi)
    }
  }

  a0=1; a1=0

  newg= newg2[!is.na(newg2)]
  newg= newg*a0 +a1

  if( length(newg) < ns/2 ) {
    return(1e6)
  }

  lam=list()
  lam= lmomgum(vec2par(c(0,1),'gum'))

  lgum=list()
  lgum=lmoms.md.park(newg, mtrim=FALSE, no.stop=TRUE)

  if(lgum$ifail == 1) {
    return(1e6)
  }

  zvec[1] = lam$lambdas[1]*a0 + a1 - lgum$lambdas[1]
  zvec[2] = lam$lambdas[2]*a0  - lgum$lambdas[2]
  zvec[3] = lam$lambdas[3]  - lgum$lambdas[3]

  z= t(zvec) %*% covinv %*% zvec

  nllh.norm =  z/2   + (3/2)*log( (2*pi) ) + lcovdet

  if(pen=='norm' | pen=="normal"){

    pk.ns = -log( pk.norm.stnary(para= a, mu= mu, std= std) )

  }else if(pen=='beta' | pen=="Beta"){

    pk.ns = -log( pk.beta.ns(para= a, lme.center=lme.center, p=p,
                             c1=c1,c2=c2)$pk.one )

  }else if(pen=="no"){
    pk.ns =0

  }else if(pen=='ms' | pen=="MS"){

    pk.ns = -log( pk.beta.stnary(para=a, p=6, q=9)$pk.one)

  }else if(pen=="park" |pen=="Park"){

    pk.ns = -log( pk.beta.stnary(para=a, p=2.5, q=2.5)$pk.one)

  }else if(pen=="cannon" | pen=="Cannon"){

    pk.ns = -log( pk.beta.stnary(para=a, p=2, q=3.3)$pk.one  )

  }else if(pen=="cd" | pen=="CD"){

    if (a[3] >= 0) {pk.ns <- 0
    }else if (a[3] > -1 & a[3] < 0) {
      pk.ns <- -log( exp(-((1/(1 + a[3])) - 1)) )
    }else if (a[3] <= -1) {pk.ns = 10^6
    }

  }

  pen.out= max(abs(xi)- 1.0, 0)

  nllh.glme = nllh.norm + pk.ns + pen.out*abs(xi)*100

  return(nllh.glme)
}

#----------------------------------------------------------------------
#' Generalized L-moments estimation for non-stationary GEV11 model
#'
#' @description
#' This function estimates parameters of the non-stationary GEV11 model
#' where mu(t) = mu0 + mu1*t and sigma(t) = exp(sigma0 + sigma1*t).
#'
#' @param xdat A numeric vector of data to be fitted.
#' @param ntry Number of attempts for parameter estimation (default 10).
#' @param ftol Tolerance for convergence (default 1e-6).
#' @param init.rob Use robust regression for initialization (default TRUE).
#' @param glme.pre Pre-estimation method: "wls" (default) or "gado".
#' @param opt.choose Selection criterion: "gof" (default, goodness-of-fit) or "nllh" (negative log-likelihood).
#' @param pen Type of penalty function: "norm", "beta" (default), "ms", "park", "cannon", "cd", or "no".
#' @param pen.choice Choice number for penalty hyperparameters (default 6 for beta).
#' @param mu Mean for normal penalty (default -0.55).
#' @param std Std for normal penalty (default 0.3).
#' @param p Shape for beta penalty (default 6).
#' @param c1 Scaling for beta penalty (default 10).
#' @param c2 Limit for beta penalty (default 5).
#'
#' @return A list containing:
#' \itemize{
#'   \item para.glme - Proposed GLME estimates (5 parameters: mu0, mu1, sigma0, sigma1, xi).
#'   \item para.lme - L-moment based estimates for non-stationary model.
#'   \item para.gado - GN16 original estimates.
#'   \item para.wls - Weighted least squares estimates (WLS).
#'   \item strup.org - WLSE by strup method.
#'   \item lme.sta - Stationary L-moment estimates.
#'   \item pen - Penalty method used.
#'   \item p_q - (for beta) p and q values.
#'   \item c1_c2 - (for beta) c1 and c2 values.
#' }
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
#'   \code{\link{nsgev}} for the pure L-moment wrapper (no penalty),
#'   \code{\link{quagev.NS}} for non-stationary quantile computation.
#'
#' @author Yonggwan Shin, Seokkap Ko, Jihong Park, Yire Shin, Jeong-Soo Park
#'
#' @examples
#' # Load example streamflow data
#' data(streamflow)
#' x <- streamflow$r1
#'
#' \donttest{
#' # Estimate non-stationary GEV11 parameters
#' result <- glme.gev11(x, ntry = 5)
#' print(result$para.glme)  # Proposed GLME estimates
#' print(result$para.lme)  # L-moment based estimates
#' }
#'
#' @export
glme.gev11 = function(xdat, ntry=10, ftol=1e-6,
                      init.rob=TRUE, glme.pre="wls",
                      opt.choose="gof", pen='beta', pen.choice=NULL,
                      mu= -0.55, std=0.3, p=6, c1=10, c2=5){

  z <- list()
  ns=length(xdat)
  year=seq(1,length(xdat))

  model='gev11'
  name_gev11_ns  =c("mu0","mu1","sigma0","sigma1","xi")
  name_gev00_sta =c("mu_sta","sigma_sta","xi_sta")

  if(is.null(pen) | !is.character(pen)) {
    stop("pen should be given as a character")}
  if(pen=="Beta" ) pen="beta"; if(pen=="CD") pen="cd"
  if(pen=="MS") pen="ms"; if(pen=="Park") pen="park"
  if(pen=="normal") pen="norm"; if(pen=="Cannon") pen="cannon"

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

  # --------------Strup WLS ---------------------------------------------
  reg.dat=data.frame( cbind(year, xdat) )

  mu.init= lm(xdat~year, reg.dat)$coefficients
  m0= mu.init[1]
  m1= mu.init[2]

  orig.para=c(m0,m1,1.0,-0.001,-0.1)

  strup = strup.glme.gev11(xdat, orig.para=orig.para, rob=FALSE)

  #---------------- GN16 ------------------------------------------------
  qlist= make.qmax.gev11(xdat, orig.para=orig.para, rob=FALSE)

  orig.para=c(m0,m1, qlist$sig0, qlist$sig1, 0)

  gado = time.m.gev11(qmax=qlist$qmax, orig.para=orig.para,
                      rob=FALSE, mdfy= FALSE)  # GN16

  #------ proposed method --------------------------------------------
  if(init.rob==TRUE){
    mu.init= lmrob(xdat~year, reg.dat)$coefficients  # robust regression
  }else{
    mu.init= lm(xdat~year, reg.dat)$coefficients #OLS reg
  }

  m0.rob= mu.init[1]
  m1.rob= mu.init[2]

  orig.para=c(m0.rob, m1.rob, 1.0,-0.001,0)

  if(glme.pre=="gado"){

    qlist= make.qmax.gev11(xdat, orig.para=orig.para, rob=TRUE)

    orig.para=c(m0.rob, m1.rob, qlist$sig0, qlist$sig1, 0)

    gado.rob = time.m.gev11(qmax=qlist$qmax, orig.para=orig.para,
                            rob=TRUE, mdfy=TRUE)

    pretheta= gado.rob$para.org

  }else if(glme.pre=="wls"){

    orig.para=c(m0.rob, m1.rob, strup$strup.final[3],
                strup$strup.final[4], 0)

    strup.glme = strup.glme.gev11(xdat, orig.para=orig.para,
                                  rob=TRUE)

    pretheta= strup.glme$strup.final
  }

  z = optim.glme.gev11(xdat, ntry=ntry, ftol=ftol,
                       pretheta=pretheta, model=model,
                       pen=pen, mu=mu, std=std,
                       p=p,c1=c1,c2=c2, choose=opt.choose)

  if(z$precis > ftol) { z$para.lme = pretheta
  message("no optim for proposed lme") }

  # --------------------------------------------------------------------
  z$para.gado   = gado$para.org            # GN16 original est
  z$strup.org   = strup$strup.para         # wlse by strup
  z$para.wls =  strup$strup.final          # specified WLSE

  z$lme.sta = pargev(lmoms(xdat,nmom=3))$para   # stationary L-ME

  names(z$para.glme)     <-name_gev11_ns
  names(z$para.lme)      <-name_gev11_ns
  names(z$para.gado)     <-name_gev11_ns
  names(z$strup.org)     <-name_gev11_ns
  names(z$para.wls)      <-name_gev11_ns
  names(z$lme.sta)       <-name_gev00_sta

  if(pen=='beta'){

    ww= pk.beta.ns(para=z$para.glme[c(1,3,5)], lme.center=z$lme.sta,
                   p=p,c1=c1,c2=c2)
    z$p_q= c(ww$p,ww$q)
    z$c1_c2=c(c1,c2)
  }

  return(z)
}

#-------------------------------------------------
#' L-moment distance function for GEV11 model
#'
#' @description Internal function that computes the L-moment distance vector
#' for the GEV11 model. The residuals are formed from the difference
#' between theoretical Gumbel L-moments and sample L-moments of
#' standardized data.
#'
#' @param a Numeric vector (mu0, sigma0, xi) to optimize.
#' @param xdat Numeric vector of data.
#' @param pretheta Full parameter vector (mu0, mu1, sigma0, sigma1, xi) from pre-estimation.
#'
#' @return Numeric vector of length 3, representing L-moment equation residuals.
#'
#' @keywords internal
obj.lme.gev11 <- function(a, xdat=xdat, pretheta=pretheta)
{
  zz=rep(100,3)

  mu0 <- a[1]
  mu1 <- pretheta[2]
  sig0 <- a[2]
  sig1 = pretheta[4]
  xi= a[3]

  ns=length(xdat)
  year=seq(1,ns)
  gum01=rep(NA, ns)
  gum.dat=rep(NA, ns)
  newg2=rep(NA,ns)

  gum.dat[1:ns]= xdat[1:ns]-(mu0 + mu1*year[1:ns])
  gum.dat[1:ns]= gum.dat[1:ns]/exp(sig0 + sig1*year[1:ns])

  gum01[1:ns]= 1-xi*gum.dat[1:ns]

  for (it in 1:ns ) {
    if( is.na(gum01[it]) ){
      newg2[it]=NA
    }else if( gum01[it] <= 0 ) {
      newg2[it]= NA
    }else if( gum01[it] > 0) {
      newg2[it]= log(gum01[it])/(-xi)
    }
  }

  a0=1; a1=0

  newg=newg2[!is.na(newg2)]
  newg= newg*a0 + a1

  if( length(newg) < ns/2 ) {
    zz[1:3]=1000
    return(zz)
  }

  lam=list()
  lam= lmomgum(vec2par(c(0,1),'gum'))

  lgum=list()
  lgum=lmoms.md.park(newg, mtrim=FALSE, no.stop=TRUE)

  if(lgum$ifail == 1) {
    zz[1:3]=1000
    return(zz)
  }

  pen.out= max(abs(xi)- 1.0, 0)

  zz[1] = lam$lambdas[1]*a0 + a1 - lgum$lambdas[1]
  zz[2] = lam$lambdas[2]*a0  - lgum$lambdas[2]
  zz[3] = lam$lambdas[3]  - lgum$lambdas[3]

  zz[3] = zz[3] + sign(zz[3])*pen.out*100
  return(zz)
}

#---------------------------------------------------
#' Multi-start optimization for GEV11 model
#'
#' @description Internal function that performs multi-start optimization
#' for the non-stationary GEV11 model. First solves L-moment equations
#' via nleqslv (Broyden method), then optionally refines with penalized
#' optimization via optim. Selects the best solution using goodness-of-fit
#' or penalized negative log-likelihood.
#'
#' @param xdat Numeric vector of data.
#' @param ntry Number of random starting points. Default is 10.
#' @param ftol Tolerance for convergence. Default is 1e-6.
#' @param pretheta Pre-estimated parameter vector (mu0, mu1, sigma0, sigma1, xi).
#' @param model Model type: "gev11" (default), "gev10", or "gev20".
#' @param pen Penalty type: "beta" (default), "norm", "ms", "park", "cannon", "cd", or "no".
#' @param mu Normal penalty mean.
#' @param std Normal penalty standard deviation.
#' @param p Beta penalty shape parameter.
#' @param c1 Beta penalty scaling parameter.
#' @param c2 Beta penalty limit parameter.
#' @param choose Selection method: "gof" (goodness-of-fit) or "nllh".
#'
#' @return A list containing:
#' \describe{
#'   \item{para.lme}{L-moment based estimates (5 parameters)}
#'   \item{precis}{Precision of the best solution}
#'   \item{para.glme}{GLME estimates (if pen != "no")}
#'   \item{nllh.glme}{Penalized negative log-likelihood (if pen != "no")}
#'   \item{pen}{Penalty type used}
#' }
#'
#' @keywords internal
optim.glme.gev11= function(xdat, ntry=10, ftol=1e-6,
                           pretheta=NULL, model="gev11",
                           pen='beta', mu=mu, std=std, p=p,
                           c1=c1, c2=c2, choose=NULL)
{
  zm=list()
  value=list()
  k=list()

  init = matrix(0, nrow=ntry, ncol=3)
  init = init.glme.gev11(xdat, ntry, pretheta)

  if(model=='gev10') npar=4
  if(model=='gev20') npar=5
  if(model=='gev11') npar=5

  precis=rep(1000, ntry)
  para.sel=matrix(NA,ntry,ncol=npar)

  tryCatch({
    for(i in 1:ntry) {
      value =  tryCatch( nleqslv( x=as.vector(init[i,1:3]),
                                  fn= obj.lme.gev11,
                                  method="Broyden",
                                  xdat=xdat, pretheta=pretheta) )
      k[[i]] <- value

      if(is(value)[1]=="try-error"){
        k[[i]]$fvec <- 10^6
        k[[i]]$termcd = 5
      }else{
        precis[i]=  mean(abs(k[[i]]$fvec) )

        if( precis[i] < ftol) {
          k[[i]]$root = value$x
          para.sel[i,1:5]=c( k[[i]]$root[1], pretheta[2],
                             k[[i]]$root[2], pretheta[4], k[[i]]$root[3])
        }
      }

      precis[is.na(precis[i])]=1000
      if( abs( k[[i]]$termcd ) > 3 ) {
        precis[i]=1000
        para.sel[i,]=NA
      }
    } #end for
  }) # trycathch

  sel.fin = sel.para_all(xdat,         #para est. SSP JKSS (2025)
                         para.sel, model, obj.fun=abs(precis))

  zm$para.lme= sel.fin$para
  zm$precis = precis[sel.fin$min.itry]

  if(pen != "no"){   # perform glme

    gev11.cov =list()

    gev11.cov =gev11.GLD(par=zm$para.lme, xdat=xdat)

    covinv = gev11.cov$covinv
    lcovdet= gev11.cov$lcovdet

    # parameter estimation using optim and glme

    gntry= min(5,ntry)
    my.nllh= rep(1e6,gntry)
    para.sel=matrix(NA,gntry,ncol=npar)

    newtheta= zm$para.lme
    init[2,1:3] = c(newtheta[1],newtheta[3],newtheta[5]-.01)

    tryCatch(
      for(i in 1:gntry){

        value=list()

        value <- try(
          optim(par=as.vector(init[i,1:3]), fn= nllh.glme.gev11,
                xdat=xdat, newtheta=newtheta, covinv=covinv,
                lcovdet=lcovdet, pen=pen,
                mu=mu, std=std, p=p, c1=c1,c2=c2)
        )

        if(is(value)[1]=="try-error"){
          k[[i]]$fvec <- 10^6
        }else{
          k[[i]] <- value
          k[[i]]$root = value$par
          k[[i]]$fvec = value$value
        }

        if( value$convergence != 0) { my.nllh[i]=10^6
        }else{
          my.nllh[i] = k[[i]]$fvec
          para.sel[i,1:5]=c( k[[i]]$root[1], newtheta[2],
                             k[[i]]$root[2], newtheta[4], k[[i]]$root[3])
        }

      } #for
    ) #tryCatch

    if(all(my.nllh==10^6)) {
      message("-- No solution was found in optim for glme --")
      return(zm)
    }

    if(choose=="nllh"){

      selc_num = which.min( my.nllh )    #my.nllh=k[[i]]$fvec

      x  <- k[[selc_num]]

      zm$para.glme = c(x$root[1],newtheta[2],x$root[2],
                       newtheta[4],x$root[3])
      zm$nllh.glme = x$fvec

    }else if(choose=="gof"){

      sel.fin = sel.para_all(xdat,    #para est. SSP JKSS (2025)
                             para.sel, model, obj.fun=my.nllh)

      zm$para.glme = sel.fin$para
      zm$nllh.glme = my.nllh[sel.fin$min.itry]
    }

  }else if(pen=="no"){
    zm$para.glme = zm$para.lme
  } #if pen

  zm$pen=pen
  return(zm)
}

#-------------------------------------------------
#' Calculate GLD covariance for GEV11 model
#'
#' @description Internal function that computes the generalized L-moment
#' distance covariance matrix for standardized residuals from the GEV11 model.
#' Uses bootstrap if the direct L-moment covariance matrix is singular.
#'
#' @param par Numeric vector of GEV11 parameters (mu0, mu1, sigma0, sigma1, xi).
#' @param xdat Numeric vector of data.
#'
#' @return A list containing:
#' \describe{
#'   \item{covinv}{3x3 inverse covariance matrix of L-moments}
#'   \item{lcovdet}{Log determinant of the covariance matrix}
#' }
#'
#' @keywords internal
gev11.GLD <- function(par=NULL, xdat=NULL)
{

  z=list()
  mu0 <- par[1]
  mu1 <- par[2]
  sig0 <- par[3]
  sig1 = par[4]
  xi= par[5]

  ns=length(xdat)
  year=seq(1,ns)
  gum01=rep(NA, ns)
  gum.dat=rep(NA, ns)
  newg2=rep(NA,ns)

  gum.dat[1:ns]= xdat[1:ns]-(mu0 + mu1*year[1:ns])
  gum.dat[1:ns]= gum.dat[1:ns]/exp(sig0 + sig1*year[1:ns])

  gum01[1:ns]= 1-xi*gum.dat[1:ns]

  for (it in 1:ns ) {
    if( is.na(gum01[it]) ){
      newg2[it]=NA

    }else if( gum01[it] <= 0 ) {
      newg2[it]= NA
    }else if( gum01[it] > 0) {
      newg2[it]= log(gum01[it])/(-xi)
    }
  }

  a0=1; a1=0

  newg= newg2[!is.na(newg2)]
  newg= newg*a0 + a1

  covinv= matrix(NA, 3,3)
  cov=lmoms.cov(newg, nmom=3)
  covinv=solve(cov)

  detc = det(cov)

  if(detc <= 0){
    BB=200          # we need Bootstrap to calculate cov ---
    sam.lmom= matrix(NA,BB,3)

    for (ib in 1:BB){
      sam.lmom[ib,1:3]=lmoms(sample(newg,size=ns,replace=TRUE),
                             nmom=3)$lambdas
    }
    cov=cov(sam.lmom)
    covinv=solve(cov)
    detc=det(cov)
  }

  z$covinv =covinv
  z$lcovdet =log(detc)

  return(z)
}

#-------------------------------------------------
#' Select best parameters based on goodness-of-fit
#'
#' @description Internal function that selects the best parameter estimate
#' from multiple candidates using energy distance-based goodness-of-fit
#' across multiple return periods.
#'
#' @param xdat Numeric vector of data.
#' @param para.sel Matrix of candidate parameter sets (ntry x npar).
#' @param model Model type: "gev11", "gev10", or "gev20".
#' @param obj.fun Numeric vector of objective function values (tie-breaker).
#'
#' @return A list containing:
#' \describe{
#'   \item{para}{Best parameter vector}
#'   \item{min.itry}{Index of the best candidate}
#'   \item{gof}{Goodness-of-fit values for all candidates}
#'   \item{obj.fun}{Input objective function values}
#' }
#'
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

#-------------------------------------------------
#' Goodness-of-fit based on exceedance counts
#'
#' @description Internal function that computes goodness-of-fit by comparing
#' observed vs expected exceedances above return level quantiles for multiple
#' return periods.
#'
#' @param xdat Numeric vector of data.
#' @param vecT Numeric vector of return periods. Default is c(5, 10, 20, 40, 80).
#' @param para Parameter vector for the model.
#' @param model Model type: "gev11", "gev10", "gev20", or "gev00".
#'
#' @return Sum of absolute relative exceedance errors (scalar).
#'
#' @keywords internal
gof.ene_all = function(xdat, vecT=c(5,10,20,40,80),
                       para=NULL, model=NULL){

  ns=length(xdat)
  nT = length(vecT)
  chi=rep(NA,nT)

  for(i in 1:nT){
    qt = quagev.NS(f=1-(1/vecT[i]), para, nsample=ns, model)
    ene = ns/vecT
    sne = sum(xdat >= qt)
    chi[i] = abs(ene[i]-sne) /ene[i]
  }
  sum(chi)
}

#--------------------------------------------------
#' Modified L-moments calculation with optional trimming
#'
#' @description Internal function that computes sample L-moments with an
#' option for trimmed L-moments (TLmoms with left trimming). Includes
#' a fail-safe mode for use in optimization loops.
#'
#' @param x Numeric vector of data.
#' @param nmom Number of L-moments to compute. Default is 5.
#' @param mtrim Logical. If TRUE, use left-trimmed L-moments (trim=5). Default is FALSE.
#' @param no.stop Logical. If TRUE, return failure indicator instead of stopping. Default is FALSE.
#' @param vecit Logical. If TRUE, return as a vector. Default is FALSE.
#'
#' @return A list (from \code{lmomco::TLmoms()}) with an additional \code{ifail}
#'   component (0 = success, 1 = failure).
#'
#' @keywords internal
lmoms.md.park = function (x, nmom = 5, mtrim=FALSE,
                          no.stop = FALSE, vecit = FALSE)
{
  z=list()
  ifail=0

  n <- length(x)
  if (nmom > n) {
    if (no.stop) {
      ifail=1
      z$ifail=ifail
      return(z)
    }else{
      stop("More L-moments requested by parameter 'nmom' than data points available in 'x'")
    }
  }
  if (length(unique(x)) == 1) {
    if (no.stop) {
      ifail=1
      z$ifail=ifail
      return(z)
    }else{stop("all values are equal--Lmoments can not be computed")
    }
  }

  if(ifail == 0) {

    if(mtrim==FALSE){
      z <- TLmoms(x, nmom = nmom)
    }else if(mtrim==TRUE){
      z <- TLmoms(x, nmom = nmom, leftrim=5)
    }
    z$source <- "lmoms"
    if (!vecit)
      z$ifail=ifail
    return(z)
    if (nmom == 1) {
      z <- z$lambdas[1]
    }
    else if (nmom == 2) {
      z <- c(z$lambdas[1], z$lambdas[2])
    }
    else {
      z <- z$lambdas[1:nmom]
    }
    attr(z, which = "trim") <- NULL
    attr(z, which = "rightrim") <- NULL
    attr(z, which = "leftrim") <- NULL
    attr(z, which = "source") <- "lmoms"

  }
  z$ifail=ifail
  return(z)
}

#-----------------------------------------------------
#' Create maximum residual series for GEV11 model
#'
#' @description Internal function that constructs a modified residual series
#' (qmax) from the data after removing estimated location trend, used for
#' the GN16 time-varying moment method.
#'
#' @param xdat Numeric vector of data.
#' @param orig.para Initial parameter vector (mu0, mu1, sigma0, sigma1, xi).
#' @param rob Logical. If TRUE, use robust regression. If FALSE, use OLS.
#'
#' @return A list containing:
#' \describe{
#'   \item{qmax}{Modified residual series (numeric vector)}
#'   \item{sig0}{Log-scale intercept estimate}
#'   \item{sig1}{Log-scale trend estimate}
#' }
#'
#' @keywords internal
make.qmax.gev11 =function(xdat=NULL, orig.para=NULL, rob=NULL)
{

  z=list()
  m0= orig.para[1]
  m1= orig.para[2]
  ns = length(xdat)
  year= seq(1,ns)

  res = xdat - (m0 + m1* year)
  mres= mean(res)
  res.pr = abs(res - mres)

  lres.pr=log(res.pr)
  sig.dat=data.frame( cbind(year, lres.pr) )

  if(rob==FALSE){
    sig.lm= lm(lres.pr~year, sig.dat)$coefficients
  }else if(rob==TRUE){
    sig.lm= lmrob(lres.pr~year, sig.dat)$coefficients
  }

  sig0 = sig.lm[1]
  sig1 = sig.lm[2]
  sigt = exp(sig0 + sig1* year)

  qmax= rep(NA, ns)

  for(i in 1:ns){
    if(sig1 >= 0){
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
  z$sig0=sig0
  z$sig1=sig1
  return(z)
}

#------------------------------------------------------
#' Initialize parameters for multi-start GEV11 optimization
#'
#' @description Internal function that generates initial parameter sets
#' (mu0, sigma0, xi) for multi-start optimization in the GEV11 model.
#' Uses stationary L-moment estimates and random perturbations.
#'
#' @param data Numeric vector of data.
#' @param ntry Number of initial parameter sets to generate. Default is 10.
#' @param pretheta Pre-estimated parameter vector (mu0, mu1, sigma0, sigma1, xi).
#'
#' @return A matrix with \code{ntry} rows and 3 columns (mu0, log(sigma), xi).
#'
#' @keywords internal
init.glme.gev11 <-function(data, ntry=10, pretheta=NULL){

  init <-matrix(0, nrow=ntry, ncol=3)
  if(abs(pretheta[5]) > 0.5) pretheta[5] = sign(pretheta[5])*0.48

  lmom_init = lmoms(data,nmom=5)
  lmom_est <- pargev(lmom_init)

  init[1,1] <- lmom_est$para[1]
  init[1,2] <- log(lmom_est$para[2])
  init[1,3] =  lmom_est$para[3]

  if( abs(lmom_est$para[3]) > 0.5) {
    init[1,3] = sign(init[1,3])*0.48 }

  maxm1=ntry-2; maxm2=ntry-3
  init[2:maxm1,1] <- init[1,1]+rnorm(n=maxm2,mean=0,sd = 20)
  init[2:maxm1,2] <- log(lmom_est$para[2])+rnorm(n=maxm2,mean=0,sd = 1)
  init[2:maxm1,3] <- runif(n=maxm2,min= -0.49, max=0.49)

  mx = mean(data)
  sx= log(sqrt(var(data)))
  init[ntry-1,1:3] = c(mx, sx, pretheta[5])
  init[ntry,1:3] =   c(pretheta[1], pretheta[3], pretheta[5]+.05)
  return(init)
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
#'   \code{\link{glme.gev}} for stationary GEV estimation.
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

#--------------------------------------------
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

#-----strup wls -----------------------------------
#' Strup WLS estimation for GEV11
#'
#' @description Internal function that performs weighted least squares (WLS)
#' estimation for the non-stationary GEV11 model using the Strup method.
#' Estimates location trend and log-scale parameters through iterative
#' regression steps.
#'
#' @param xdat Numeric vector of data.
#' @param orig.para Initial parameter vector (mu0, mu1, sigma0, sigma1, xi).
#' @param rob Logical. If TRUE, use robust regression. If FALSE, use OLS.
#'
#' @return A list containing:
#' \describe{
#'   \item{strup.sta}{Stationary L-moment estimates of standardized residuals}
#'   \item{strup.para}{Raw Strup parameter estimates (5 parameters)}
#'   \item{strup.final}{Final adjusted parameter estimates (5 parameters)}
#' }
#'
#' @keywords internal
strup.glme.gev11 =function(xdat, orig.para=NULL, rob=NULL){

  w=list()
  ns=length(xdat)
  year=seq(1,ns)

  m0=orig.para[1]
  m1=orig.para[2]
  res=xdat -(m0+m1*year)

  stand= wls.gev11(xdat, res, rob=rob)

  new.para= c(stand$m, stand$sig, 0)
  ares= stand$res
  sigt =  exp(new.para[3]+new.para[4]*year)

  w$strup.sta = pargev(lmoms(ares))$para

  strup.para = c(new.para[1:4], w$strup.sta[3])

  #------to specify final parameter values -----------------------------
  yt=rep(0,ns+1)
  mu_st=w$strup.sta[1]
  year2=seq(0,ns)

  for (ka in 1:(ns+1) ) {
    ti= ka-1
    yt[ka]=  strup.para[1] + strup.para[2] * ti
    yt[ka]= yt[ka] + mu_st * exp(strup.para[3]+ strup.para[4]* ti)
  }

  nh=round((ns/2))
  yt[nh-1] = yt[nh-1] + 0.02;   yt[nh-2] = yt[nh-2] - 0.02
  yt[nh+1] = yt[nh+1] - 0.01;   yt[nh+2] = yt[nh+2] + 0.01

  reg.dat=data.frame( cbind(year2, yt) )

  if(rob==FALSE){
    mu.init= lm(yt ~ year2, reg.dat)$coefficients
  }else if(rob==TRUE){
    mu.init= lmrob(yt ~ year2, reg.dat)$coefficients
  }

  sigmaf_0 = strup.para[3] + log(w$strup.sta[2])  #w$strup.sta[2] = sig_st
  sigmaf_1 = strup.para[4]
  xif = strup.para[5]

  w$strup.para = strup.para
  w$strup.final= c(mu.init, sigmaf_0, sigmaf_1, xif)
  return(w)
}

#----------------------------------------------------------
#' Weighted least squares core estimation for GEV11
#'
#' @description Internal function that performs the core WLS steps:
#' (Step 3) estimate log-scale trend from absolute residuals,
#' (Step 4) estimate location trend from weighted data,
#' (Step 5) compute standardized residuals.
#'
#' @param xdat Numeric vector of data.
#' @param res Numeric vector of residuals from initial location regression.
#' @param rob Logical. If TRUE, use robust regression. If FALSE, use OLS.
#'
#' @return A list containing:
#' \describe{
#'   \item{sig}{Log-scale regression coefficients c(sigma0, sigma1)}
#'   \item{m}{Location regression coefficients c(mu0, mu1)}
#'   \item{res}{Standardized residuals}
#' }
#'
#' @keywords internal
wls.gev11 = function(xdat, res=NULL, rob=NULL){

  z=list()
  ns=length(res)
  year=seq(1,ns)

  lres.pr=log(abs(res))
  sig.dat=data.frame( cbind(year, lres.pr) )

  if(rob==FALSE){
    z$sig= lm(lres.pr~year, sig.dat)$coefficients             # Step 3
  }else if(rob==TRUE){
    z$sig= lmrob(lres.pr~year, sig.dat)$coefficients
  }

  sigt =  exp(z$sig[1] + z$sig[2]* year)

  res.n = xdat/sigt
  ytran0 = rep(1,ns)/sigt
  ytran1 = year/sigt

  new.data= data.frame( cbind(ytran0, ytran1, res.n) )

  if(rob==FALSE){
    z$m= lm(res.n~0+ytran0+ytran1, new.data)$coefficients      # Step 4
  }else if(rob==TRUE){
    z$m= lmrob(res.n~0+ytran0+ytran1, new.data)$coefficients
  }

  z$res = res.n -(z$m[1]*ytran0+z$m[2]*ytran1)               # Step 5
  return(z)
}

#------------------------------------------------------
#' Time-varying moment estimation (GN16 method)
#'
#' @description Internal function that implements the GN16 time-varying
#' moment method for non-stationary GEV11 estimation. Computes time-varying
#' location from the shape parameter and scale function.
#'
#' @param qmax Modified residual series from \code{make.qmax.gev11()}.
#' @param orig.para Initial parameter vector (mu0, mu1, sigma0, sigma1, xi).
#' @param rob Logical. If TRUE, use robust regression. Default is FALSE.
#' @param mdfy Logical. If TRUE, apply modified GN16 with updated alpha0. Default is TRUE.
#'
#' @return A list containing:
#' \describe{
#'   \item{mu_t}{Time-varying location values}
#'   \item{para.org}{Estimated parameters (mu0, mu1, sigma0, sigma1, xi)}
#'   \item{mu_t.up}{Updated time-varying location (if mdfy=TRUE)}
#' }
#'
#' @keywords internal
time.m.gev11 = function(qmax=NULL, orig.para=NULL, rob=FALSE,
                        mdfy=TRUE){

  z=list()
  para.gado=rep(NA,5)
  ns=length(qmax)
  year=seq(1,ns)

  m0=orig.para[1]
  m1=orig.para[2]
  sig0=orig.para[3]
  sig1=orig.para[4]

  lmom_q = lmoms(qmax)
  q.sta = pargev(lmom_q)$para
  xi = q.sta[3]

  cd = sqrt( (xi^2) /( gamma(1+2*xi) - gamma(1+xi)^2 ) )
  alpha_t = exp(sig0 +sig1*year) * cd
  z$mu_t=  - (1-gamma(1+ xi) )*alpha_t/xi  + m0+ m1*year

  #--------------------------------------------------------
  nh=round((ns/2))
  mu.gado=z$mu_t
  mu.gado[nh-1] = mu.gado[nh-1] + 0.02
  mu.gado[nh+1] = mu.gado[nh+1] - 0.02

  mu.data= data.frame( cbind(year, mu.gado) )
  if(rob==TRUE){
    loc.gado =lmrob(mu.gado~year, mu.data)$coefficients
  }else if(rob==FALSE){
    loc.gado =lm(mu.gado~year, mu.data)$coefficients
  }

  alpha0 = log(cd) + sig0
  alpha1 = sig1

  if(mdfy != TRUE){
    z$para.org= c(loc.gado, alpha0, alpha1, xi )

  }else if(mdfy==TRUE){
    # modify GN16 -----------------------------------------
    alpha0_up = log(q.sta[2])- sig1*nh
    z$para.org = c(loc.gado, alpha0_up, alpha1, xi )

    alpha_t = exp(alpha0_up +sig1*year)
    z$mu_t.up=  - (1-gamma(1+ xi) )*alpha_t/xi  + m0+ m1*year
    # -----------------------------------------------------
  }
  return(z)
}
