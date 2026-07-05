# Model Averaging for GEV High Quantile Estimation
#
# Reference: Shin, Y., Shin, Y., & Park, J. S. (2026). Model averaging with
# mixed criteria for estimating high quantiles of extreme values: Application
# to heavy rainfall. Stochastic Environmental Research and Risk Assessment,
# 40(2), 47. https://doi.org/10.1007/s00477-025-03167-x


#' Model Averaging for GEV High Quantile Estimation
#'
#' @description
#' This function estimates high quantiles of the Generalized Extreme Value (GEV)
#' distribution using model averaging with mixed criteria. It combines Maximum
#' Likelihood Estimation (MLE) and L-moment Estimation (LME) to construct
#' candidate submodels and assign weights effectively.
#'
#' @param data A numeric vector of data to be fitted (e.g., annual maxima).
#' @param quant The probabilities corresponding to high quantiles to be estimated.
#'   Default is c(0.98, 0.99, 0.995).
#' @param weight The weighting method name. Options are:
#'   \itemize{
#'     \item 'like', 'like0', 'like1' (default): Likelihood-based weights (AIC)
#'     \item 'gLd', 'gLd0', 'gLd1', 'gLd2': Generalized L-moment distance weights
#'     \item 'med', 'med1', 'med2': Median-based weights
#'     \item 'cvt': Conventional AIC weights
#'   }
#'   Variants with numbers indicate left trimming level (0, 1, or 2).
#' @param numk The number of candidate submodels K. Default is 12.
#' @param B The number of bootstrap samples. Default is 200.
#' @param varcom Logical. Whether to compute variance of quantile estimates.
#'   Default is TRUE.
#' @param trim The number of left trimming for L-moments. Usually 0 (default), 1, or 2.
#' @param fig Logical. Whether to produce diagnostic plots. Default is FALSE.
#' @param bma Logical. Whether to use Bayesian Model Averaging. Default is FALSE.
#' @param pen Penalty type for BMA prior: 'norm' (normal, default) or 'beta'.
#' @param CD Logical. Whether to compute Coles-Dixon penalized MLE. Default is FALSE.
#' @param remle Logical. Whether to compute restricted MLE. Default is FALSE.
#'
#' @details
#' The model averaging approach works as follows:
#' \enumerate{
#'   \item MLE and LME of GEV parameters are computed
#'   \item K candidate shape parameters (xi) are selected from profile likelihood CI
#'   \item For each candidate xi, MLE with fixed xi is computed
#'   \item Weights are assigned based on the selected method
#'   \item Final quantile estimates are weighted averages across submodels
#' }
#'
#' The weighting schemes include:
#' \itemize{
#'   \item 'like': AIC-based weights using likelihood with fixed xi
#'   \item 'gLd': Weights based on generalized L-moment distance
#'   \item 'med': Weights based on median and L-moment distance
#'   \item 'cvt': Conventional AIC weights
#' }
#'
#' When bma=TRUE, Bayesian model averaging is applied with prior specified by pen.
#'
#' @return A list containing:
#' \itemize{
#'   \item mle.hosking - MLE estimates in Hosking style (mu, sigma, xi)
#'   \item qua.mle - Quantile estimates from MLE
#'   \item mle.cov3 - Covariance matrix of MLE (3x3)
#'   \item qua.se.mle.delta - Standard errors of MLE quantiles (delta method)
#'   \item lme - L-moment estimates (mu, sigma, xi)
#'   \item lme.cov3 - Covariance matrix of LME (bootstrap)
#'   \item qua.lme - Quantile estimates from LME
#'   \item qua.se.lme.boots - Standard errors of LME quantiles (bootstrap)
#'   \item qua.ma - Model-averaged quantile estimates
#'   \item w.ma - Weights used for model averaging
#'   \item fixw.se.ma - Asymptotic SE under fixed weights
#'   \item ranw.se.ma - Asymptotic SE under random weights
#'   \item surr - Surrogate model parameters (mu, sigma, xi)
#'   \item pick_xi - Selected xi values for K submodels
#'   \item qua.bma - (if bma=TRUE) BMA quantile estimates
#'   \item w.bma - (if bma=TRUE) BMA weights
#'   \item mle.CD - (if CD=TRUE) Coles-Dixon penalized MLE
#'   \item qua.CD - (if CD=TRUE) Quantile estimates from CD-MLE
#'   \item remle1 - (if remle=TRUE) Restricted MLE (first constraint)
#'   \item qua.remle1 - (if remle=TRUE) Quantile estimates from remle1
#'   \item remle2 - (if remle=TRUE) Restricted MLE (second constraint)
#'   \item qua.remle2 - (if remle=TRUE) Quantile estimates from remle2
#'   \item quant - The quantile probabilities used
#' }
#'
#' @references
#' Shin, Y., Shin, Y., & Park, J. S. (2026). Model averaging with mixed criteria
#' for estimating high quantiles of extreme values: Application to heavy rainfall.
#' \emph{Stochastic Environmental Research and Risk Assessment}, 40(2), 47.
#' \doi{10.1007/s00477-025-03167-x}
#'
#' @seealso \code{\link{glme.gev}} for stationary GLME estimation,
#'   \code{\link{magev.ksensplot}} for K sensitivity analysis,
#'   \code{\link{magev.qqplot}} for Q-Q diagnostic plots,
#'   \code{\link{magev.rlplot}} for return level plots.
#'
#' @author Yonggwan Shin, Seokkap Ko, Jihong Park, Yire Shin, Jeong-Soo Park
#'
#' @examples
#' # Load example data
#' data(haenam)
#' x <- haenam$X1
#'
#' # Basic usage with likelihood weights
#' result <- ma.gev(x, quant = c(0.95, 0.99), weight = 'like1', B = 100)
#' print(result$qua.ma)     # Model-averaged quantiles
#' print(result$qua.mle)    # MLE quantiles for comparison
#' print(result$qua.lme)    # LME quantiles for comparison
#'
#' # Using generalized L-moment distance weights
#' result2 <- ma.gev(x, quant = c(0.95, 0.99), weight = 'gLd', B = 100)
#' print(result2$w.ma)      # Model weights
#'
#' @export
ma.gev <- function(data = NULL, quant = c(0.98, 0.99, 0.995),
                   weight = 'like1', numk = 12, B = 200, varcom = TRUE,
                   trim = 0, fig = FALSE, bma = FALSE, pen = "norm",
                   CD = FALSE, remle = FALSE) {

  zx <- list()
  SMALL <- 1e-05
  numq <- length(quant)
  zx$weight <- weight
  start <- "mle"
  zx$start <- start
  zx$quant <- quant

  # Internal defaults
  pick <- 0.95
  pertr <- 0.85
  surr <- TRUE
  numom <- 3
  order <- 3
  boot.lme <- TRUE
  cov.lme <- NULL
  type <- 'full'
  cov.type <- 'lambda'
  fix.kpar <- FALSE
  pick.xi <- NULL

  if (fix.kpar == FALSE) {
    pick.xi <- NULL
  } else if (fix.kpar == TRUE) {
    numk <- length(pick.xi)
  }
  if (surr == FALSE) fig <- FALSE
  if (fig == TRUE) {
    surr <- TRUE
    boot.lme <- TRUE
    varcom <- TRUE
  }

  # Normalize weight names
  if (weight == 'gld') weight <- 'gLd'
  if (cov.type == 'r') cov.type <- 'ratio'
  if (cov.type == 'l') cov.type <- 'lambda'
  if (weight == 'like0') {
    weight <- 'like'
    trim <- 0
  }
  if (weight == 'like1') {
    weight <- 'like'
    trim <- 1
  }
  if (weight == 'gld0' | weight == 'gLd0') {
    weight <- 'gLd'
    trim <- 0
  }
  if (weight == 'gld1' | weight == 'gLd1') {
    weight <- 'gLd'
    trim <- 1
  }
  if (weight == 'gld2' | weight == 'gLd2') {
    weight <- 'gLd'
    trim <- 2
  }
  if (weight == 'med1') {
    weight <- 'med'
    trim <- 1
  }
  if (weight == 'med2') {
    weight <- 'med'
    trim <- 2
  }

  zx$data <- data
  org.numk <- numk

  para3 <- matrix(NA, nrow = numk, ncol = 3)
  wtgd <- rep(NA, numk)
  zp <- matrix(NA, numq, numk)
  zpf <- rep(NA, numq)

  nsample <- length(data)
  delta <- list()

  # ------- MLE and SE by delta method -------------------
  delta <- gev.rl.delta(data = data, ntry = 10, quant = quant)

  zx$mle.hosking <- delta$mle
  zx$qua.mle <- delta$qua.mle
  zx$mle.cov3 <- delta$cov

  if (varcom == TRUE & !is.null(delta$cov)) {
    qua.se.mle <- rep(NA, numq)
    for (i in 1:numq) {
      qua.se.mle[i] <- delta.gev(gevf3 = delta, quant = quant[i], d3yes = TRUE)$v3
    }
    zx$qua.se.mle.delta <- sqrt(qua.se.mle)
  }

  # ----- LME and SE by bootstrap --------------------
  if (is.null(cov.lme)) {
    if (start == 'lme' | weight == 'med' | start == 'mix' |
        weight == 'lcv' | weight == 'opt') {
      boot.lme <- TRUE
    }
  }

  hosking <- list()
  hosking <- lme.boots(data = data, B = B, quant = quant,
                           boot = boot.lme, trim = trim)

  if (boot.lme == TRUE) {
    zx$lme.cov3 <- hosking$cov.par
  } else if (boot.lme == FALSE) {
    zx$lme.cov3 <- cov.lme$cov.par
  }

  zx$lme <- hosking$lme
  zx$qua.lme <- hosking$qua.lme
  zx$leftrim <- trim
  hosking$mle <- delta$mle

  xi.hat <- (delta$mle[3] + hosking$lme[3]) / 2

  # ---- CD penalized MLE and restricted MLE ------------------------------------------
  if (CD == TRUE) {
    CD.mle <- rep(NA, 3)
    CD.mle <- mle.gev.CD(xdat = data, ntry = 5)$mle  # hosking style xi

    if (any(is.na(CD.mle))) {
      zx$qua.CD <- NA
      zx$mle.CD <- NA
      CD.mle <- NULL
    } else {
      zx$qua.CD <- quagev(quant, vec2par(CD.mle, 'gev'))
      zx$mle.CD <- CD.mle
    }
  }

  if (remle == TRUE) {
    if (CD != TRUE) CD.mle <- hosking$lme
    rem.try <- list()
    rem.try <- remle.gev(xdat = data, ntry = 5, rest = 'mean', quant = quant,
                         CD.mle = CD.mle, mle = delta$mle,
                         second = TRUE, w.mpse = FALSE)

    zx$remle1 <- rem.try$para.remle1
    zx$qua.remle1 <- rem.try$qua.remle1
    zx$remle2 <- rem.try$para.remle2
    zx$qua.remle2 <- rem.try$qua.remle2
  }
  # --------------------------------------------------------------

  hosking$start <- start
  hosking$numk <- numk
  hosking$quant <- quant
  hosking$weight <- weight
  hosking$B <- B
  hosking$trim <- trim
  hosking$data <- data
  hosking$boot.lme <- boot.lme

  # --- xi_k picking for candidate submodels ---------------------------
  kpar <- rep(NA, numk)

  if (fix.kpar == FALSE) {
    kpar.list <- list()

    kpar.list <- cand.xi(data, hosking = hosking, mle = zx$mle.hosking,
                                   pick0 = pick, nint = 256, start = start,
                                   numk = numk, figure = fig, bma = bma, pen = pen)

    kpar <- kpar.list$kpar
    hosking$start <- kpar.list$start
  } else if (fix.kpar == TRUE) {
    kpar <- pick.xi
    hosking$start <- start
  }
  hosking$kpar <- kpar
  hosking$cov.type <- cov.type

  # ------- Zp and weights computing -----------------------------------
  run <- 1

  xqa <- c(0.5, .65, .8, .85, .9, .925, .95, .965, .98,
           .985, .99, .9925, .995, .9965, .998, .999)
  mywt <- list()
  hosking$pertr <- pertr
  bmaw <- rep(0, numk)
  wtgd <- rep(NA, numk)

  mywt <- weight.com(data, numk = numk, hosking = hosking,
                         kpar = kpar, numom = numom,
                         xqa = xqa, varcom = varcom, boot.lme = boot.lme,
                         cov.lme = cov.lme, surr = surr,
                         type = type, trim = trim, cov.type = cov.type,
                         bma = bma, pen = pen)

  para3 <- mywt$prob.call$mle3
  npara3 <- na.omit(para3)
  numk <- nrow(npara3)
  if (mywt$numk != numk) message("Note: numk adjusted after weight computation")
  notid <- which(!is.na(para3[, 1]))
  mywt$data <- data
  para3 <- npara3
  mywt$prob.call$mle3 <- npara3
  kpar <- kpar[notid]

  zp <- mywt$zp[, notid, drop = FALSE]
  xzp <- mywt$xzp[, notid, drop = FALSE]
  wtgd <- mywt$wtgd[notid]
  mywt$wtgd <- wtgd
  bmaw <- mywt$bmaw[notid]

  if (bma != TRUE) {
    newk <- new.kpar2(wtgd = wtgd, numk = numk, kpar = kpar)
  } else if (bma == TRUE) {
    newk <- new.kpar2(wtgd = bmaw, numk = numk, kpar = kpar)
  }

  numk <- newk$numk
  aw <- newk$aw

  while (aw >= 1) {
    run <- run + 1
    kpar <- newk$kpar2
    hosking$kpar <- newk$kpar2
    hosking$numk <- numk

    mywt <- list()
    bmaw <- rep(0, numk)
    wtgd <- rep(NA, numk)
    para3 <- matrix(NA, nrow = numk, ncol = 3)
    zp <- matrix(NA, numq, numk)

    mywt <- weight.com(data, numk = numk, hosking = hosking,
                           kpar = newk$kpar2, numom = numom,
                           xqa = xqa, varcom = varcom, boot.lme = boot.lme,
                           cov.lme = cov.lme, surr = surr,
                           type = type, trim = trim, cov.type = cov.type,
                           bma = bma, pen = pen)

    para3 <- mywt$prob.call$mle3
    npara3 <- na.omit(para3)
    numk <- nrow(npara3)
    if (mywt$numk != numk) message("Note: numk adjusted in iteration")
    notid <- which(!is.na(para3[, 1]))
    mywt$data <- data
    para3 <- npara3
    mywt$prob.call$mle3 <- npara3
    kpar <- mywt$kpar
    kpar <- kpar[notid]

    zp <- mywt$zp[, notid, drop = FALSE]
    xzp <- mywt$xzp[, notid, drop = FALSE]
    wtgd <- mywt$wtgd[notid]
    mywt$wtgd <- wtgd
    bmaw <- mywt$bmaw[notid]
    bmaw2 <- if (!is.null(mywt$bmaw2)) mywt$bmaw2[notid] else NULL

    if (run >= 5) break
    mxid <- which.max(wtgd)
    bmxid <- which.max(bmaw)

    tre1 <- 0.1
    tre8 <- 0.5

    if (bma != TRUE) {
      if (wtgd[1] > tre1 | wtgd[numk] > tre1 | wtgd[mxid] > tre8) {
        newk <- new.kpar2(wtgd = wtgd, numk = numk, kpar = kpar)
        numk <- newk$numk
        aw <- newk$aw
      } else {
        aw <- 0
        break
      }
    }

    if (bma == TRUE) {
      if (bmaw[1] > tre1 | bmaw[numk] > tre1 | bmaw[bmxid] > tre8) {
        newk <- new.kpar2(wtgd = bmaw, numk = numk, kpar = kpar)
        numk <- newk$numk
        aw <- newk$aw
      } else {
        aw <- 0
        break
      }
    }
  }

  idp3 <- which(is.na(para3[, 1]))
  if (length(idp3) > 0) {
    para3[idp3, ] <- 0
    wtgd[idp3] <- 0
    bmaw[idp3] <- 0
  }

  # Handle NA in zp - different handling for numq==1 vs numq>1
  if (numq > 1) {
    id <- which(is.na(t(zp)[, 1]))
    if (length(id) > 0) {
      zp[, id] <- 0
      wtgd[id] <- 0
      bmaw[id] <- 0
    }
  } else if (numq == 1) {
    id <- which(is.na(zp))
    if (length(id) > 0) {
      zp[id] <- 0
      wtgd[id] <- 0
      bmaw[id] <- 0
    }
  }

  if (all(wtgd == 0)) {
    message("All return levels are NA, returning LME quantiles")
    zx$qua.ma <- zx$qua.lme
    zx$qua.bma <- zx$qua.lme
    class(zx) <- "magev"
    return(zx)
  }

  bmaw[which(is.na(bmaw))] <- 0.0
  wtgd[which(is.na(wtgd))] <- 0.0

  zpf <- rep(NA, numq)
  zpf.bma <- rep(NA, numq)
  if (numq > 1) {
    for (iq in 1:numq) {
      zpf[iq] <- t(wtgd) %*% t(zp)[1:numk, iq]
      if (bma == TRUE) zpf.bma[iq] <- t(bmaw) %*% t(zp)[1:numk, iq]
    }
  } else if (numq == 1) {
    zpf[1] <- sum(wtgd * as.vector(zp))
    if (bma == TRUE) zpf.bma[1] <- sum(bmaw * as.vector(zp))
  }

  # ---- calculating asymptotic SE by delta method -----------------
  if (varcom == TRUE) {
    covint <- cov.interp(numk, para3, cov2 = mywt$prob.call$cov2)
    avar <- asymp.var(mywt, covint, qqq = quant, order)
  }

  # ----- finding surrogate for MA using GEV -----------
  if (surr == TRUE) {
    zpf.surr <- t(wtgd) %*% t(xzp)
    para.ma <- t(wtgd) %*% para3
    surro <- surrogate(zpf.surr, xqa, init.surr = para.ma)
    zx$surr <- list()
    zx$surr$par <- as.vector(surro$par)
  }

  # ----- Predictive posterior variance ----------------
  if (varcom == TRUE & bma == TRUE) {
    zpdiff <- matrix(NA, numq, numk)
    loc.var <- matrix(NA, numq, numk)
    for (ip in 1:numk) {
      if (numq > 1) {
        zpdiff[1:numq, ip] <- (zp[1:numq, ip] - zpf.bma[1:numq])^2
      } else if (numq == 1) {
        zpdiff[1:numq, ip] <- (zp[ip] - zpf.bma[1:numq])^2
      }
      loc.var[1:numq, ip] <- avar$MatC[ip, ip, 1:numq]
    }
    msm <- t(bmaw) %*% t(zpdiff)
    mse <- t(bmaw) %*% t(loc.var)
  }

  if (varcom == TRUE) {
    zx$fixw.se.ma <- avar$fin.se.MA.qua
    zx$ranw.se.ma <- avar$adj.se.MA.qua
    if (bma == TRUE) {
      zx$pred.se.bma <- as.vector(sqrt(msm + mse))
      zx$bma.se.between <- sqrt(msm)
      zx$bma.se.within <- sqrt(mse)
    }
  }

  if (boot.lme == TRUE) {
    zx$qua.se.lme.boots <- hosking$qua.lme.se
  } else {
    zx$qua.se.lme.boots <- cov.lme$qua.lme.se
  }

  # ------- finalization ---------------------------
  if (bma != TRUE) {
    zx$zp <- zp
    zx$qua.ma <- zpf
    zx$w.ma <- wtgd
    zx$pick_xi <- kpar
    zx$run.numk <- numk
  }

  if (bma == TRUE) {
    zx$qua.bma <- zpf.bma
    zx$pen <- pen
    zx$w.bma <- bmaw

    set <- set.prior(pen = pen, numk = numk, xi_lme = hosking$lme[3],
                     kpar = kpar, weight = weight)
    if (pen == "norm") zx$prior_mu_std <- set$prior_mu_std
    if (pen == "beta") zx$p_q_beta <- set$p_q_beta

    zx$pick_xi <- kpar
    zx$run.numk <- numk
  }

  zx$run_kpar <- run
  zx$original.numk <- org.numk

  class(zx) <- "magev"
  return(zx)
}
