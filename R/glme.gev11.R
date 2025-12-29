# R code for Non-Stationary GEV11 model using L-moment
# GEV11: mu(t) = mu0 + mu1*t, sigma(t) = exp(sigma0 + sigma1*t)

#-------------------------------------------------------------------------
#' Normal preference function for non-stationary GEV
#'
#' @param para A vector of GEV parameters.
#' @param mu Mean for normal distribution.
#' @param std Standard deviation for normal distribution.
#' @return Preference function value.
#' @author Jeong-Soo Park
#' @keywords internal
pk.norm.ns = function(para=NULL, mu=NULL, std=NULL){
  Brone = 1 + dnorm(para[3], mean= mu, sd=std)
  return(Brone)
}

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
#' @author Jeong-Soo Park
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

  Bef <- function(x) { ((-al+x)^(p-1)) * ((bl-x)^(q-1)) }
  Be  <- integrate(Bef, lower=al, upper=bl)[1]$value

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
#' @author Jeong-Soo Park
#' @keywords internal
gev.glme.m0s0_11 <- function(a, xdat=xdat, newtheta=newtheta,
                             covinv=covinv, lcovdet=lcovdet,
                             pen=pen, mu=mu, std=std, p=p, c1=c1, c2=c2)
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
    pk.ns = -log( pk.norm.ns(para= a, mu= mu, std= std) )
  }else if(pen=='beta' | pen=="Beta"){
    pk.ns = -log( pk.beta.ns(para= a, lme.center=lme.center, p=p,
                             c1=c1,c2=c2)$pk.one )
  }else if(pen=="no"){
    pk.ns =0
  }else if(pen=='ms' | pen=="MS"){
    pk.ns = -log( MS_pk(para=a, p=6, q=9))
  }else if(pen=="park" |pen=="Park"){
    pk.ns = -log( MS_pk(para=a, p=2.5, q=2.5))
  }else if(pen=="cannon" | pen=="Cannon"){
    pk.ns = -log( MS_pk(para=a, p=2, q=3.3)  )
  }else if(pen=="cd" | pen=="CD"){
    if (a[3] >= 0) {pk.ns <- 0
    }else if (a[3] > -1 & a[3] < 0) {
      pk.ns <- -log( exp(-((1/(1 + a[3])) - 1)) )
    }else if (a[3] <= -1) {pk.ns = 10^6
    }
  }

  nllh.glme = nllh.norm + pk.ns

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
#' @param pen Type of penalty function: "norm", "beta" (default), "ms", "park", "cannon", "cd", or "no".
#' @param pen.choice Choice number for penalty hyperparameters.
#' @param mu Mean for normal penalty (default -0.55).
#' @param std Std for normal penalty (default 0.3).
#' @param p Shape for beta penalty (default 6).
#' @param c1 Scaling for beta penalty (default 10).
#' @param c2 Limit for beta penalty (default 5).
#'
#' @return A list containing:
#' \itemize{
#'   \item para.glme - Proposed GLME estimates (5 parameters: mu0, mu1, sigma0, sigma1, xi).
#'   \item para.jkss - L-moment based estimates for non-stationary model.
#'   \item para.gado - GN16 original estimates.
#'   \item strup.sta - Stationary WLSE.
#'   \item strup.org - WLSE by strup.
#'   \item strup.final - Specified WLSE.
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
#' @author Jeong-Soo Park
#'
#' @examples
#' # Load example streamflow data
#' data(streamflow)
#' x <- streamflow$r1
#'
#' # Estimate non-stationary GEV11 parameters
#' \donttest{
#' result <- glme.gev11(x, ntry = 5)
#' print(result$para.glme)  # Proposed GLME estimates
#' print(result$para.jkss)  # L-moment based estimates
#' }
#'
#' @export
glme.gev11 = function(xdat, ntry=10, ftol=1e-6, init.rob=TRUE,
                      pen='beta', pen.choice=NULL, mu=-0.55, std=0.3,
                      p=6, c1=10, c2=5){

  z <- list()
  ns=length(xdat)
  year=seq(1,length(xdat))

  model='gev11'
  name_gev11_ns  =c("mu0","mu1","sigma0","sigma1","xi")
  name_gev00_sta =c("mu_sta","sigma_sta","xi_sta")

  if(pen=="Beta" ) pen="beta"; if(pen=="CD") pen="cd"
  if(pen=="MS") pen="ms"; if(pen=="Park") pen="park"
  if(pen=="normal") pen="norm"; if(pen=="Cannon") pen="cannon"

  if(pen=='beta' & !is.null(pen.choice)){
    pc1c2= matrix(c(6,6,6,2,2,2,10,20,30,10,20,30,5,7,9,5,7,9),
                  6,3, byrow=FALSE)
    p= pc1c2[pen.choice,1]
    c1=pc1c2[pen.choice,2]
    c2=pc1c2[pen.choice,3]
  }
  if(pen=='norm' & !is.null(pen.choice)){
    mustd= matrix(c(-.5,.2,-.5,.1,-.6,.2,-.6,.1),4,2, byrow=TRUE)
    mu= mustd[pen.choice,1]
    std=mustd[pen.choice,2]
  }

  # --------------Strup WLS ---------------------------------------------
  reg.dat=data.frame( cbind(year, xdat) )

  mu.init= lmrob(xdat~year, reg.dat)$coefficients      # robust regression
  m0= mu.init[1]
  m1= mu.init[2]

  orig.para=c(m0,m1,1.0,-0.001,-0.1)

  strup = strup_11_glme_const(xdat, orig.para=orig.para, const=FALSE)

  #---------------- GN16 ------------------------------------------------
  qlist= make.qmax_11(xdat, orig.para=orig.para, rob=FALSE)

  orig.para=c(m0,m1, qlist$sig0, qlist$sig1, 0)

  gado = calc_time_m_11(qmax=qlist$qmax, orig.para=orig.para)  # GN16

  #------ proposed method --------------------------------------------
  if(init.rob==TRUE){
    mu.init= lmrob(xdat~year, reg.dat)$coefficients  # robust regression
  }else{
    mu.init= lm(xdat~year, reg.dat)$coefficients #OLS reg
  }

  m0.rob= mu.init[1]
  m1.rob= mu.init[2]

  orig.para=c(m0.rob, m1.rob, 1.0,-0.001,0)

  qlist= make.qmax_11(xdat, orig.para=orig.para, rob=TRUE)

  orig.para=c(m0.rob, m1.rob, qlist$sig0, qlist$sig1, 0)

  gado.rob = calc_time_m_11(qmax=qlist$qmax, orig.para=orig.para)

  z = multi.m0s0_11(xdat, ntry=ntry, ftol=ftol,
                    pretheta=gado.rob$para.org, model=model,
                    pen=pen, mu=mu, std=std,
                    p=p,c1=c1,c2=c2)

  if(z$precis > ftol) { z$para.jkss = gado.rob$para.org
  cat("no optim for proposed","\n") }

  # --------------------------------------------------------------------
  z$para.gado   = gado$para.org            # GN16 orginal est
  z$strup.sta   = strup$strup.sta          # stationary wlse
  z$strup.org   = strup$strup.para         # wlse by strup
  z$strup.final = strup$strup.final        # specified WLSE

  z$lme.sta = pargev(lmoms(xdat,nmom=5))$para   # stationary L-ME

  names(z$para.glme)     <-name_gev11_ns
  names(z$para.jkss)     <-name_gev11_ns
  names(z$para.gado)     <-name_gev11_ns
  names(z$strup.org)     <-name_gev11_ns
  names(z$strup.final)   <-name_gev11_ns
  names(z$strup.sta)     <-name_gev00_sta
  names(z$lme.sta)       <-name_gev00_sta

  if(pen=='beta'){
    ww=pk.beta.ns(para=z$para.glme[c(1,3,5)], lme.center=z$lme.sta,
                  p=p, c1=c1,c2=c2)
    z$p_q= c(ww$p,ww$q)
    z$c1_c2=c(c1,c2)
  }

  return(z)
}

#-------------------------------------------------
#' L-moment distance function for GEV11 model
#' @keywords internal
gev.Ldist.m0s0_11 <- function(a, xdat=xdat, pretheta=pretheta)
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

  pen= max(abs(xi)- 1.0, 0)

  zz[1] = lam$lambdas[1]*a0 + a1 - lgum$lambdas[1]
  zz[2] = lam$lambdas[2]*a0  - lgum$lambdas[2]
  zz[3] = lam$lambdas[3]  - lgum$lambdas[3]

  zz[3] = zz[3] + sign(zz[3])*pen
  return(zz)
}

#---------------------------------------------------
#' Multi-start optimization for GEV11 model
#' @keywords internal
multi.m0s0_11= function(xdat, ntry=10, ftol=1e-6,
                        pretheta=pretheta, model=model,
                        pen='beta', mu=mu, p=p, std=std, c1=c1, c2=c2)
{
  zm=list()
  value=list()
  k=list()

  init = matrix(0, nrow=ntry, ncol=3)
  init = ginit.m0s0(xdat, ntry, pretheta)

  if(model=='gev10') npar=4
  if(model=='gev20') npar=5
  if(model=='gev11') npar=5

  precis=rep(1000, ntry)
  para.sel=matrix(NA,ntry+1,ncol=npar)

  tryCatch({
    for(i in 1:ntry) {
      value =  tryCatch( nleqslv( x=as.vector(init[i,1:3]),
                                  fn= gev.Ldist.m0s0_11,
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
                             k[[i]]$root[2], pretheta[4],  k[[i]]$root[3])
        }
      }

      precis[is.na(precis[i])]=1000
      if( abs( k[[i]]$termcd ) > 3 ) {
        precis[i]=1000
        para.sel[i,]=NA
      }
    }
  })

  zm$para.jkss =sel.para_all(xdat, para.sel, model)$para  # L-moment estimates
  zm$precis =precis[which.min(precis)]

  if(pen != "no"){   #  perform glme

    gntry= 3
    gev11.cov =list()
    isol = 0

    gev11.cov =gev11.GLD(par=zm$para.jkss, xdat=xdat)

    covinv = gev11.cov$covinv
    lcovdet= gev11.cov$lcovdet

    my.nllh=rep(1e6,gntry)

    newtheta= zm$para.jkss
    init[2,1:3] = c(newtheta[1],newtheta[3],newtheta[5]-.01)

    tryCatch(
      for(i in 1:gntry){

        value=list()

        value <- try(
          optim(par=as.vector(init[i,1:3]), fn= gev.glme.m0s0_11,
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

        if( value$convergence != 0) {my.nllh[i]=10^6
        }else{
          isol=isol+1
          my.nllh[i] = k[[i]]$fvec
        }

      }
    )

    if(isol==0) {
      cat("-- No solution was found in nleqslv or optim --","\n")
      return(zm)
    }

    selc_num = which.min( my.nllh )

    x  <- k[[selc_num]]

    zm$nllh.glme = k[[selc_num]]$fvec
    zm$para.glme = c(x$root[1],newtheta[2],x$root[2],newtheta[4],x$root[3])

  }

  # When pen="no", set para.glme same as para.jkss
  if(pen == "no"){
    zm$para.glme = zm$para.jkss
  }

  zm$pen=pen
  return(zm)
}

#-------------------------------------------------
#' Calculate GLD covariance for GEV11 model
#' @keywords internal
gev11.GLD <- function(par=NULL, xdat=xdat)
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
    BB=200
    sam.lmom= matrix(NA,BB,3)

    for (ib in 1:BB){
      sam.lmom[ib,1:3]=lmoms(sample(newg,size=ns,replace=T), nmom=3)$lambdas
    }
    cov=cov(sam.lmom)
    covinv=solve(cov)
    detc=det(cov)
  }

  lcovdet =log(detc)
  z$covinv =covinv
  z$lcovdet =lcovdet

  return(z)
}

#-------------------------------------------------
#' Select best parameters based on GOF
#' @keywords internal
sel.para_all =function(xdat, para.sel, model=NULL){

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
    z$para =upara.sel[length(gof),]
  }else{
    z$para =upara.sel[which.min(gof),]
  }

  z$gof=gof
  return(z)
}

#-------------------------------------------------
#' Goodness-of-fit function
#' @keywords internal
gof.ene_all = function(xdat, vecT=c(5,10,20,40,80), para, model=NULL){

  ns=length(xdat)
  nT = length(vecT)
  year=seq(1,ns)
  chi=rep(NA,nT)

  for(i in 1:nT){
    Tp =  vecT[i]
    qt = qns.gev_all(Tp, para, year, model)
    ene = ns/Tp
    sne = sum(xdat >= qt)
    chi[i] = abs(ene-sne) /ene
  }
  chi2 = sum(chi)
  return(chi2)
}

#--------------------------------------------------
#' Modified L-moments calculation
#' @keywords internal
lmoms.md.park =
  function (x, nmom = 5, mtrim=FALSE, no.stop = FALSE, vecit = FALSE)
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
      }else if(mtrim==T){
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
#' Create qmax for GEV11 model
#' @keywords internal
make.qmax_11 =function(xdat=NULL, orig.para=NULL, rob=NULL)
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
#' Initialize parameters for multi-start
#' @keywords internal
ginit.m0s0 <-function(xdat, ntry=ntry, pretheta){

  init <-matrix(0, nrow=ntry, ncol=3)
  if(abs(pretheta[5]) > 0.5) pretheta[5] = sign(pretheta[5])*0.48

  lmom_init = lmoms(xdat,nmom=5)
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

  mx = mean(xdat)
  sx= log(sqrt(var(xdat)))
  init[ntry-1,1:3] = c(mx, sx, pretheta[5])
  init[ntry,1:3] =   c(pretheta[1], pretheta[3], pretheta[5]+.05)
  return(init)
}

#----------------------------------------------------
#' Quantile function for non-stationary GEV
#' @keywords internal
qns.gev_all= function(Tp=NULL, para, year, model=NULL){

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
#' Set parameters based on model type
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
  }

  z$mu=c(mu0,mu1,mu2)
  z$sig=c(sig0,sig1)
  z$xi=xi
  return(z)
}

#---------------------------------------------------------
#' Quantile function for GEV110 model
#' @keywords internal
qns.gev110= function(Tp, para, year){

  nsample=length(year)
  ns=nsample
  zpT=rep(NA, nsample)

  mu0=para[1]
  mu1=para[2]
  sigma0=para[3]
  sigma1=para[4]
  xi=para[5]

  zpc= (1- ( -log(1-(1/Tp) ) )^xi ) /xi
  zpT[1:ns]= mu0+mu1*year[1:ns] + zpc*exp(sigma0 + sigma1*year[1:ns])

  return(zpT)
}

#-----strup wls -----------------------------------
#' Strup WLS with robust regression
#' @keywords internal
strup_11_glme_const_lmrob =function(xdat, orig.para=orig.para,
                                    const=NULL){

  w=list()
  ns=length(xdat)
  year=seq(1,ns)

  m0=orig.para[1]
  m1=orig.para[2]
  res=xdat -(m0+m1*year)

  stand= wls.park_11_lmrob(xdat,res)

  new.para= c(stand$m, stand$sig, 0)
  ares= stand$res
  sigt =  exp(new.para[3]+new.para[4]*year)

  w$strup.sta = pargev(lmoms(ares))$para

  if(const==TRUE){
    if( abs(w$strup.sta[3]) > 0.5 ){
      xifix= sign(w$strup.sta[3])*0.4999
      work= pargev.kfix(lmoms(ares, nmom=3), kfix= xifix)
      w$strup.sta = work$para
    } }

  strup.para = c(new.para[1:4], w$strup.sta[3])

  yt=rep(0,ns+1)
  mu_st=w$strup.sta[1]
  year2=seq(0,ns)

  for (ka in 1:(ns+1) ) {
    t= ka-1
    yt[ka]=  strup.para[1] + strup.para[2] * t
    yt[ka]= yt[ka] + mu_st * exp(strup.para[3]+ strup.para[4]* t)
  }

  nh=round((ns/2))
  yt[nh-1] = yt[nh-1] + 0.03;   yt[nh-2] = yt[nh-2] - 0.02
  yt[nh+1] = yt[nh+1] - 0.03;   yt[nh+2] = yt[nh+2] + 0.02

  reg.dat=data.frame( cbind(year2, yt) )

  mu.init= lm(yt ~ year2, reg.dat)$coefficients

  sigmaf_0 = strup.para[3] + log(w$strup.sta[2])
  sigmaf_1 = strup.para[4]
  xif = strup.para[5]

  w$strup.para = strup.para
  w$strup.final= c(mu.init, sigmaf_0, sigmaf_1, xif)

  return(w)
}

#-----strup wls -----------------------------------
#' Strup WLS with ordinary regression
#' @keywords internal
strup_11_glme_const =function(xdat, orig.para=orig.para, const=NULL){

  w=list()
  ns=length(xdat)
  year=seq(1,ns)

  m0=orig.para[1]
  m1=orig.para[2]
  res=xdat -(m0+m1*year)

  stand= wls.park_11(xdat,res)

  new.para= c(stand$m, stand$sig, 0)
  ares= stand$res
  sigt =  exp(new.para[3]+new.para[4]*year)

  w$strup.sta = pargev(lmoms(ares))$para

  if(const==TRUE){
    if( abs(w$strup.sta[3]) > 0.5 ){
      xifix= sign(w$strup.sta[3])*0.4999
      work= pargev.kfix(lmoms(ares, nmom=3), kfix= xifix)
      w$strup.sta = work$para
    } }

  strup.para = c(new.para[1:4], w$strup.sta[3])

  yt=rep(0,ns+1)
  mu_st=w$strup.sta[1]
  year2=seq(0,ns)

  for (ka in 1:(ns+1) ) {
    ti= ka-1
    yt[ka]=  strup.para[1] + strup.para[2] * ti
    yt[ka]= yt[ka] + mu_st * exp(strup.para[3]+ strup.para[4]* ti)
  }

  nh=round((ns/2))
  yt[nh-1] = yt[nh-1] + 0.03;   yt[nh-2] = yt[nh-2] - 0.02
  yt[nh+1] = yt[nh+1] - 0.03;   yt[nh+2] = yt[nh+2] + 0.02

  reg.dat=data.frame( cbind(year2, yt) )

  mu.init= lm(yt ~ year2, reg.dat)$coefficients

  sigmaf_0 = strup.para[3] + log(w$strup.sta[2])
  sigmaf_1 = strup.para[4]
  xif = strup.para[5]

  w$strup.para = strup.para
  w$strup.final= c(mu.init, sigmaf_0, sigmaf_1, xif)
  return(w)
}

#----------------------------------------------------------
#' WLS with robust regression
#' @keywords internal
wls.park_11_lmrob = function(xdat,res ){

  z=list()
  ns=length(res)
  year=seq(1,ns)

  lres.pr=log(abs(res))
  sig.dat=data.frame( cbind(year, lres.pr) )
  z$sig= lmrob(lres.pr~year, sig.dat)$coefficients

  sigt =  exp(z$sig[1] + z$sig[2]* year)

  res.n = xdat/sigt
  ytran0 = rep(1,ns)/sigt
  ytran1 = year/sigt

  new.data= data.frame( cbind(ytran0, ytran1, res.n) )

  z$m= lmrob(res.n~0+ytran0+ytran1, new.data)$coefficients

  z$res = res.n -(z$m[1]*ytran0+z$m[2]*ytran1)
  return(z)
}

#------------------------------------------------------
#' WLS with ordinary regression
#' @keywords internal
wls.park_11 = function(xdat,res ){

  z=list()
  ns=length(res)
  year=seq(1,ns)

  lres.pr=log(abs(res))
  sig.dat=data.frame( cbind(year, lres.pr) )
  z$sig= lm(lres.pr~year, sig.dat)$coefficients

  sigt =  exp(z$sig[1] + z$sig[2]* year)

  res.n = xdat/sigt
  ytran0 = rep(1,ns)/sigt
  ytran1 = year/sigt

  new.data= data.frame( cbind(ytran0, ytran1, res.n) )

  z$m= lm(res.n~0+ytran0+ytran1, new.data)$coefficients

  z$res = res.n -(z$m[1]*ytran0+z$m[2]*ytran1)
  return(z)
}

#------------------------------------------------------
#' Time-varying moment estimation
#' @keywords internal
calc_time_m_11 = function(qmax=NULL, orig.para=NULL){

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
  if( xi <= -0.5) xi = -0.4999

  cd = sqrt( (xi^2) /( gamma(1+2*xi) - gamma(1+xi)^2 ) )
  alpha_t = exp(sig0 +sig1*year) * cd
  z$mu_t=  - (1-gamma(1+ xi) )*alpha_t/xi  + m0+ m1*year

  nh=round((ns/2))
  mu.gado=z$mu_t
  mu.gado[nh-1] = mu.gado[nh-1] + 0.02
  mu.gado[nh+1] = mu.gado[nh+1] - 0.02

  mu.data= data.frame( cbind(year, mu.gado) )
  loc.gado =lm(mu.gado~year, mu.data)$coefficients

  alpha0 = log(cd) + sig0
  alpha1 = sig1

  z$para.org= c(loc.gado, alpha0, alpha1, xi )

  return(z)
}
