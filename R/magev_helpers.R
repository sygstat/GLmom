# Internal helper functions for MAGEV (Model Averaging GEV)
# These functions are not exported and are used internally by ma.gev()
#
# Reference: Shin, Y., Shin, Y., & Park, J.-S. (2025). Model averaging with
# mixed criteria for estimating high quantiles of extreme values.
# arXiv preprint arXiv:2505.21417.


#' @keywords internal
gev.max <- function(xdat, ntry = 5) {
  z <- list()
  k <- list()
  n <- ntry

  nsample <- length(xdat)
  z$nsample <- nsample

  init <- matrix(0, nrow = ntry, ncol = 3)
  init <- ginit.max(xdat, ntry)
  upsig <- sqrt(var(xdat)) * 5
  upmu <- abs(mean(xdat)) * 3

  gev.lik.max <- function(a) {
    mu <- a[1]
    sc <- a[2]
    xi <- a[3]

    y <- (xdat - mu) / sc
    y <- 1 - xi * y  # Hosking style

    for (i in 1:nsample) {
      y[i] <- max(0, y[i], na.rm = TRUE)
    }

    if (any(y <= 0) || any(sc <= 0))
      return(10^6)

    if (abs(xi) >= 10^(-5)) {
      ooxi <- 1 / xi
    } else {
      ooxi <- sign(xi) * 10^5
    }

    zz <- nsample * (log(sc)) + sum(exp(ooxi * log(y))) + sum(log(y) * (1 - (ooxi)))
    return(zz)
  }

  tryCatch(
    for (i in 1:nrow(init)) {
      value <- try(solnp(init[i, ], fun = gev.lik.max,
                         LB = c(-upmu, 0, -1), UB = c(upmu, upsig, 1),
                         control = list(trace = 0, outer.iter = 10,
                                        delta = 1.e-6, inner.iter = 200, tol = 1.e-5)))

      if (is(value)[1] == "try-error") {
        k[[i]] <- list(value = 10^6)
      } else {
        k[[i]] <- value
      }
    }
  )

  optim_value <- data.frame(num = 1:n, value = sapply(k, function(x) x$value[which.min(x$value)]))
  optim_table1 <- optim_value[order(optim_value$value), ]
  selc_num <- optim_table1[1, "num"]

  x <- k[[selc_num]]

  z$conv <- x$convergence
  z$nllh <- x$value[which.min(x$value)]
  z$mle <- x$par

  return(z)
}


#' @keywords internal
ginit.max <- function(data = NULL, ntry = NULL) {
  if (ntry < 2) ntry <- 2
  n <- ntry
  init <- matrix(NA, ntry, 3)

  lmom_init <- lmoms(data, nmom = 3)
  lmom_est <- pargev(lmom_init)

  init[1, 1] <- lmom_est$para[1]
  init[1, 2] <- lmom_est$para[2]
  init[1, 3] <- lmom_est$para[3]

  maxm1 <- ntry
  maxm2 <- maxm1 - 1
  init[2:maxm1, 1] <- init[1, 1] + rnorm(n = maxm2, mean = 0, sd = 5)
  init[2:maxm1, 2] <- init[1, 2] + rnorm(n = maxm2, mean = 3, sd = 3)
  init[2:maxm1, 3] <- runif(n = maxm2, min = -0.45, max = 0.4)
  init[2:maxm1, 2] <- pmax(0.1, init[2:maxm1, 2])

  return(init)
}


#' @keywords internal
pargev.xifix.ma <- function(lmom, xifix = 0.1, checklmom = TRUE, ...) {
  para <- rep(NA, 3)
  G <- xifix
  para[3] <- G

  if (length(lmom$lambdas[1]) == 0) {
    lmom <- lmorph(lmom)
  }

  SMALL <- 1e-05
  if (abs(G) < SMALL) {
    para[1:2] <- pargum(lmom)$para
    return(list(type = "gev", para = para, source = "pargev"))
  }

  GAM <- exp(lgamma(1 + G))
  para[2] <- lmom$lambdas[2] * G / (GAM * (1 - 2^(-G)))
  para[1] <- lmom$lambdas[1] - para[2] * (1 - GAM) / G

  return(list(type = "gev", para = para, source = "pargev"))
}


#' @keywords internal
gev.rl.delta_new <- function(data, ntry = 5, quant) {
  zg <- list()
  numq <- length(quant)

  kmle.h <- gev.max(data, ntry = ntry)

  kmle.coles <- gev.fit(data, show = FALSE)

  if (kmle.h$nllh < kmle.coles$nllh) {
    if (abs(kmle.h$mle[3]) >= 0.49) kmle.h$mle[3] <- sign(kmle.h$mle[3]) * 0.49
    better <- tryCatch(gev.fit(data, show = FALSE, muinit = kmle.h$mle[1],
                               siginit = kmle.h$mle[2], shinit = -kmle.h$mle[3]))

    better.mle <- better$mle
    better.mle[3] <- -better$mle[3]
  } else {
    better <- kmle.coles
    better.mle <- kmle.coles$mle
    better.mle[3] <- -kmle.coles$mle[3]
  }

  if (abs(better.mle[3]) >= 1) {
    better.mle[3] <- sign(better.mle[3]) * 0.98
  }

  savem <- vec2par(better.mle, 'gev')
  zg$nllh <- better$nllh
  zg$mle <- better.mle
  zg$qua.mle <- quagev(quant[1:numq], savem)
  zg$data <- data
  zg$cov <- better$cov
  zg$quant <- quant

  return(zg)
}


#' @keywords internal
ginit.xifix <- function(data, ntry) {
  initx <- matrix(0, nrow = ntry, ncol = 2)

  lmom_init <- lmoms(data, nmom = 5)
  lmom_est <- pargum(lmom_init, checklmom = TRUE)

  initx[1, 1] <- lmom_est$para[1]
  initx[1, 2] <- lmom_est$para[2]

  maxm1 <- ntry
  maxm2 <- maxm1 - 1
  initx[2:maxm1, 1] <- initx[1, 1] + rnorm(n = maxm2, mean = 0, sd = 7)
  initx[2:maxm1, 2] <- initx[1, 2] + rnorm(n = maxm2, mean = 2, sd = 2)
  initx[2:maxm1, 2] <- pmax(0.1, initx[2:maxm1, 2])

  return(initx)
}


#' @keywords internal
gev.xifix.sing <- function(xdat = NULL, xifix = -0.1, ntry = 5, varcom = NULL) {
  zx <- list()
  kx <- list()
  value <- list()

  init.xi <- matrix(0, nrow = ntry, ncol = 2)
  init.xi <- ginit.xifix(xdat, ntry)
  xi <- xifix

  upsig <- diff(quantile(xdat, probs = c(0.25, 0.75))) * 5
  upmu <- abs(median(xdat)) * 3

  tryCatch(
    for (itry in 1:nrow(init.xi)) {
      value <- try(optim(init.xi[itry, ], fn = gev.xilik2,
                         lower = c(-upmu, 0), upper = c(upmu, upsig),
                         method = "L-BFGS-B",
                         xifix = xifix, xdat = xdat))

      if (is(value)[1] == "try-error") {
        kx[[itry]] <- list(value = 10^6)
      } else {
        kx[[itry]] <- value
      }
    }
  )

  optim_value <- data.frame(num = 1:ntry, value = sapply(kx, function(x) x$value[which.min(x$value)]))
  optim_table1 <- optim_value[order(optim_value$value), ]
  selc_num <- optim_table1[1, "num"]

  x <- kx[[selc_num]]

  zx$conv <- x$convergence
  zx$nllh <- x$value[which.min(x$value)]

  if (zx$conv != 0) {
    zx$mle <- NA
  } else {
    zx$mle <- x$par
    zx$mle[3] <- xifix
    if (varcom == TRUE) {
      Hess <- PrescottW(par2 = x$par, xifix = xifix, nsam = length(xdat))
      if (is.na(Hess[1, 1])) {
        zx$cov <- matrix(NA, 2, 2)
      } else if (det(Hess) <= 0) {
        zx$cov <- matrix(NA, 2, 2)
      } else {
        zx$cov <- solve(Hess)
      }
    }
  }

  class(zx) <- "gev.xifix"
  return(zx)
}


#' @keywords internal
PrescottW <- function(par2 = NULL, xifix = NULL, nsam = NULL) {
  sig <- par2[2]
  xi <- xifix
  Hess <- matrix(NA, 2, 2)

  if (xi < 0.5) {
    pw <- (1 - xi)^2 * gamma(1 - 2 * xi)
    Hess[1, 1] <- nsam * pw / (sig^2)
    Hess[2, 2] <- nsam * (1 - 2 * gamma(2 - xi) + pw) / (sig^2 * xi^2)
    Hess[1, 2] <- nsam * (pw - gamma(2 - xi)) / (sig^2 * xi)
    Hess[2, 1] <- Hess[1, 2]
  }

  return(Hess)
}


#' @keywords internal
cons.MatC <- function(mywt, cov22, quant) {
  numk <- length(mywt$wtgd)
  numq <- length(quant)
  MatC <- array(NA, c(numk, numk, numq))
  cov2i <- matrix(NA, 2, 2)
  cov2j <- matrix(NA, 2, 2)

  para3 <- mywt$prob.call$mle3

  for (i in 1:numk) {
    mle3i <- para3[i, ]
    cov2i <- as.matrix(cov22[, , i])

    for (j in i:numk) {
      mle3j <- para3[j, ]
      cov2j <- as.matrix(cov22[, , j])

      if (any(is.na(cov2i)) | any(is.na(cov2j))) {
        MatC[i, j, 1:numq] <- NA
      } else {
        MatC[i, j, 1:numq] <- delta.gev(gevf3 = NULL, mle3i = mle3i, cov2i = cov2i,
                                        mle3j = mle3j, cov2j = cov2j,
                                        quant = quant, d3yes = FALSE)$covij[1:numq]
        MatC[j, i, 1:numq] <- MatC[i, j, 1:numq]
      }
    }
  }

  wtgd <- mywt$wtgd
  bmaw <- mywt$bmaw
  fin.se <- rep(NA, numq)
  fin.se.bma <- rep(NA, numq)

  for (iq in 1:numq) {
    fin.se[iq] <- sqrt(t(wtgd) %*% MatC[, , iq] %*% wtgd)
    fin.se.bma[iq] <- sqrt(t(bmaw) %*% MatC[, , iq] %*% bmaw)
  }

  zz <- list()
  zz$MatC <- MatC
  zz$fin.se <- fin.se
  zz$fin.se.bma <- fin.se.bma
  zz$numk <- numk
  zz$numq <- numq
  zz$wtgd <- wtgd
  zz$bmaw <- bmaw

  return(zz)
}


#' @keywords internal
gev.xilik2 <- function(a, xifix = xifix, xdat = xdat) {
  gev.xilik(a, xifix = xifix, xdat = xdat)
}


#' @keywords internal
gev.xilik <- function(a, xifix = NULL, xdat = NULL) {
  mu <- a[1]
  sc <- a[2]
  xi <- xifix

  if (sc <= 0) return(10^6)

  nsam <- length(xdat)
  y <- rep(0, nsam)
  y <- (xdat - mu) / sc
  y <- 1 - xi * y  # Hosking style

  for (i in 1:nsam) {
    y[i] <- max(0, y[i], na.rm = TRUE)
  }

  if (any(y <= 0)) return(10^6)

  if (abs(xi) >= 10^(-5)) {
    ooxi <- 1 / xi
  } else {
    ooxi <- sign(xi) * 10^5
  }

  b2 <- sum(log(y)) * (1 - ooxi)
  b3 <- sum(exp(ooxi * log(y)))

  nllh <- nsam * (log(sc)) + b2 + b3
  return(nllh)
}


#' @keywords internal
lme.boots.new <- function(data, B = NULL, quant, boot = TRUE, trim = NULL) {
  z <- list()
  numq <- length(quant)
  sam.bootL <- list()
  lme.boot <- matrix(NA, nrow = B, ncol = 3)
  ratios <- matrix(NA, nrow = B, ncol = 3)
  lme.seboot <- matrix(NA, nrow = B, ncol = numq)
  med <- rep(NA, B)
  lambdas <- matrix(NA, nrow = B, ncol = 3)

  nsam <- length(data)
  re <- sort(data)
  lmt <- lmoms(data, nmom = 3)
  klmet <- pargev(lmt, checkmom = TRUE)
  z$lme <- klmet$para[1:3]
  savet <- vec2par(klmet$para[1:3], 'gev')
  z$qua.lme <- quagev(quant[1:numq], savet)
  z$quant <- quant

  if (boot == TRUE) {
    for (ib in 1:B) {
      sam.bootL[[ib]] <- sample(data, nsam, replace = TRUE)
      lm.boot <- lmoms(sam.bootL[[ib]], nmom = 3)

      lme.boot[ib, 1:3] <- pargev(lm.boot, checklmom = FALSE)$para[1:3]
      save <- vec2par(lme.boot[ib, 1:3], 'gev')
      lme.seboot[ib, 1:numq] <- quagev(quant[1:numq], save)

      lambdas[ib, 1:3] <- lm.boot$lambdas[1:3]
      ratios[ib, 1] <- lm.boot$lambdas[1]
      ratios[ib, 2:3] <- lm.boot$ratios[2:3]
      med[ib] <- median(sam.bootL[[ib]])
    }

    lme.bb <- ratios
    lme.bb[1:B, 1] <- med[1:B]
    z$cov.med.rat <- cov(lme.bb[, 1:3])

    lme.bb <- lambdas
    lme.bb[1:B, 1] <- med[1:B]
    z$cov.med.lam <- cov(lme.bb[, 1:3])

    z$cov.lambda <- cov(lambdas[, 1:3])
    z$B <- B
    z$lme.boot <- lme.boot
    z$cov.rat <- cov(ratios[, 1:3])
    z$cov.par <- cov(lme.boot[, 1:3])

    qua.lme.se <- rep(NA, numq)
    for (qi in 1:numq) {
      qua.lme.se[qi] <- sqrt(var(lme.seboot[1:B, qi], na.rm = TRUE))
    }
    z$qua.lme.se <- qua.lme.se
    z$qua.lme.boot <- lme.seboot
  }

  return(z)
}


#' @keywords internal
rcd <- function(xi) {
  minxi <- -0.2
  nxi <- length(xi)
  pfxi <- rep(NA, nxi)

  for (i in 1:nxi) {
    if (is.na(xi[i]) | is.null(xi[i])) {
      pfxi[i] <- 0
    } else {
      if (-1.0 < xi[i] & xi[i] <= minxi) {
        b1 <- -1 / (1 + minxi)
        b0 <- 1 - b1 * minxi
        pfxi[i] <- b0 + b1 * xi[i]
      }
    }
    if (xi[i] <= -1.0 | xi[i] > minxi) pfxi[i] <- 1
  }

  return(pfxi)
}


#' @keywords internal
emp.prior <- function(mle3, mle, lme) {
  z <- list()
  xi <- mle3[, 3]
  nxi <- length(xi)
  pfxi <- rep(NA, nxi)

  mnx <- min(mle[3], lme[3])
  mxx <- max(mle[3], lme[3])
  u <- -0.6
  u0 <- -0.0

  if (mxx <= u) {
    pfxi <- 1
    return(pfxi)
  }
  if (mnx <= u) {
    mnx <- u + 0.05
    b <- 0
  } else {
    b <- 1 / (mnx - u - u0)
  }
  b2 <- 1 / (mnx - u - u0)

  for (i in 1:nxi) {
    if (is.na(xi[i]) | is.null(xi[i])) {
      pfxi[i] <- 0
    } else {
      if (u + u0 < xi[i] & xi[i] < mnx) {
        pfxi[i] <- (xi[i] - u - u0) * b
      } else if (mnx <= xi[i] & xi[i] <= mxx) {
        pfxi[i] <- 1
      } else if (mxx < xi[i] & xi[i] < mxx + mnx - u) {
        pfxi[i] <- 1 + (mxx - xi[i]) * b2
      } else if (xi[i] <= u + u0 | xi[i] >= mxx + mnx - u) {
        pfxi[i] <- 0
      }
    }
  }

  return(pfxi)
}


#' @keywords internal
movave <- function(order, numk, numq, zp, wt = NULL) {
  za <- list()
  q <- order
  ez <- matrix(NA, numq, numk)
  ew <- rep(NA, numk)

  if (numk >= (q + 1)) {
    for (ip in 1:numq) {
      ez[ip, ] <- rollmean(zp[ip, 1:numk], q, fill = "extend", align = "center")
    }
  } else if (numk <= q) {
    for (ip in 1:numq) {
      for (ik in 1:numk) {
        ez[ip, ik] <- mean(zp[ip, 1:numk])
      }
    }
  }

  if (numk >= (q + 1)) {
    ew <- rollmean(wt, q, fill = "extend", align = "center")
  } else if (numk <= q) {
    for (ik in 1:numk) {
      ew[ik] <- mean(wt[1:numk])
    }
  }

  za$ez <- ez
  za$ew <- ew
  return(za)
}


#' @keywords internal
delta.gev <- function(gevf3 = NULL, mle3i = NULL, cov2i = NULL,
                      mle3j = NULL, cov2j = NULL,
                      quant = NULL, d3yes = FALSE) {
  z <- list()
  numq <- length(quant)
  p <- 1 - quant
  yp <- -log(1 - p)

  if (d3yes == TRUE) {
    c3 <- gevf3$cov
    mle <- gevf3$mle
    xi <- mle[3]
    sig <- mle[2]
    d12 <- cbind(rep(1, numq), (1 - yp^(xi)) / xi)
    d3 <- sig * (1 - yp^(xi)) / (xi^2) + sig * (yp^(xi)) * log(yp) / xi
    delta3 <- cbind(d12, d3)

    for (iq in 1:numq) {
      z$v3[iq] <- delta3[iq, ] %*% c3 %*% delta3[iq, ]
    }
    return(z)
  }

  xi2i <- mle3i[3]
  xi2j <- mle3j[3]

  qqc <- c(seq(0.1, 0.9, 0.1), 0.95, 0.98, 0.99, 0.995, 0.998, 0.999)
  zp.cori <- quagev(qqc, vec2par(mle3i, 'gev'))
  zp.corj <- quagev(qqc, vec2par(mle3j, 'gev'))
  corr.rt <- cor(zp.cori, zp.corj)

  covij <- rep(NA, length(quant))
  vi <- rep(NA, numq)
  vj <- rep(NA, numq)
  deltai <- matrix(NA, 2, numq)
  deltaj <- matrix(NA, 2, numq)

  for (id in 1:length(quant)) {
    deltai[1:2, id] <- c(1, (1 - yp[id]^(xi2i)) / xi2i)
    deltaj[1:2, id] <- c(1, (1 - yp[id]^(xi2j)) / xi2j)

    vi[id] <- t(deltai[, id]) %*% cov2i %*% deltai[, id]
    vj[id] <- t(deltaj[, id]) %*% cov2j %*% deltaj[, id]

    covij[id] <- corr.rt * sqrt(vi[id] * vj[id])
  }

  z$covij <- covij
  return(z)
}


#' @keywords internal
surrogate <- function(zpf.surr = NULL, xqa, init.surr) {
  press.fun <- function(par) {
    theta <- par
    if (theta[2] <= 0) return(10^10)
    if (theta[3] <= -1.0 | theta[3] >= 1) return(10^10)
    mzp <- quagev(xqa, vec2par(theta, 'gev'))
    fun <- sum((zpf.surr - mzp)^2)
    return(fun)
  }

  surrx <- list()
  surrx <- optim(par = init.surr, fn = press.fun)

  surrx$zp.surrmodel <- quagev(xqa, vec2par(surrx$par, 'gev'))
  surrx$zp.MA <- zpf.surr
  return(surrx)
}


#' @keywords internal
cov.interp.new <- function(numk, para3, cov2 = NULL) {
  cov22 <- array(NA, c(2, 2, numk))

  for (ik in 1:numk) {
    cov22[, , ik] <- cov2[[ik]]
  }

  id <- which(!is.na(cov22[1, 1, ]))
  id.na <- which(is.na(cov22[1, 1, ]))
  if (length(id.na) == numk) stop("all cov22 are NA")

  if (length(id) < numk) {
    xsp <- para3[id, 3]
    y1 <- cov22[1, 1, id]
    y2 <- cov22[2, 2, id]
    y12 <- cov22[1, 2, id]

    covint <- matrix(NA, 3, numk)
    covint[1, id.na] <- spline(x = xsp, y1, method = "natural", xout = para3[id.na, 3])$y
    covint[2, id.na] <- spline(x = xsp, y2, method = "natural", xout = para3[id.na, 3])$y
    covint[3, id.na] <- spline(x = xsp, y12, method = "natural", xout = para3[id.na, 3])$y

    cov22[1, 1, id.na] <- covint[1, id.na]
    cov22[2, 2, id.na] <- covint[2, id.na]
    cov22[1, 2, id.na] <- covint[3, id.na]
    cov22[2, 1, id.na] <- covint[3, id.na]
  }

  return(cov22)
}


#' @keywords internal
asymp.var <- function(mywt, covint, qqq = NULL, order = 2) {
  zx <- list()
  consC <- list()
  wtgd <- mywt$wtgd
  bmaw <- mywt$bmaw
  numq <- length(qqq)
  numk <- length(wtgd)
  zp <- mywt$zp

  consC <- cons.MatC(mywt, covint, quant = qqq)

  zx$fin.se.MA.qua <- consC$fin.se
  zx$fin.se.bma.qua <- consC$fin.se.bma

  ezvar <- rep(NA, numq)
  ewvar <- rep(NA, numq)
  ewvar.bma <- rep(NA, numq)
  trCD <- rep(0, numq)
  trCD.bma <- rep(0, numq)

  D <- cov.dir(wtgd)
  Dbma <- cov.dir(bmaw)

  move <- movave(order, numk, numq, zp, wt = wtgd)
  ez <- move$ez
  ew <- move$ew

  moveb <- movave(order, numk, numq, zp, wt = bmaw)
  ewbma <- moveb$ew

  for (ip in 1:numq) {
    ezvar[ip] <- t(ez)[, ip] %*% D %*% ez[ip, ]
    ewvar[ip] <- t(ew) %*% consC$MatC[, , ip] %*% ew
    ewvar.bma[ip] <- t(ewbma) %*% consC$MatC[, , ip] %*% ewbma

    trCD[ip] <- sum(diag(consC$MatC[, , ip] %*% D)) + ezvar[ip]
    trCD.bma[ip] <- sum(diag(consC$MatC[, , ip] %*% Dbma)) + ezvar[ip]
  }

  zx$adj.se.MA.qua <- sqrt(ewvar + trCD)
  zx$adj.se.bma.qua <- sqrt(ewvar.bma + trCD.bma)
  zx$MatC <- consC$MatC

  return(zx)
}


#' @keywords internal
cov.dir <- function(wtgd) {
  numk <- length(wtgd)
  D <- matrix(NA, numk, numk)

  for (i in 1:numk) D[i, i] <- wtgd[i] * (1 - wtgd[i]) / 2

  for (i in 1:(numk - 1)) {
    for (j in (i + 1):numk) {
      D[i, j] <- -wtgd[i] * wtgd[j] / 2
      D[j, i] <- D[i, j]
    }
  }

  return(D)
}


# ============================================================
# Coles-Dixon penalized MLE and Restricted MLE functions
# ============================================================

#' Coles-Dixon Penalized MLE for GEV
#'
#' Computes maximum penalized likelihood estimates for GEV parameters
#' using the Coles and Dixon (1999) prior penalty on the shape parameter.
#'
#' @param xdat Numeric vector of data.
#' @param ntry Number of random starting points for optimization. Default is 10.
#'
#' @return A list containing:
#'   \item{mle}{MLE estimates (mu, sigma, xi) in Hosking style}
#'   \item{nllh}{Negative log-likelihood at the optimum}
#'   \item{conv}{Convergence status (0 = success)}
#'   \item{nsample}{Sample size}
#'
#' @references
#' Coles, S., & Dixon, M. (1999). Likelihood-based inference for extreme
#' value models. Extremes, 2(1), 5-23.
#'
#' @keywords internal
gev1.CD <- function(xdat, ntry = 10) {
  z <- list()
  k <- list()
  n <- ntry

  nsample <- length(xdat)
  z$nsample <- nsample

  init <- matrix(0, nrow = ntry, ncol = 3)
  init <- ginit.max(xdat, ntry)

  # Likelihood with Coles-Dixon penalty
  gev1.lik.CD.h <- function(a) {
    mu <- a[1]
    sc <- a[2]
    xi <- a[3]

    y <- (xdat - mu) / sc
    y <- 1 - xi * y  # Hosking style

    for (i in 1:nsample) {
      y[i] <- max(0, y[i], na.rm = TRUE)
    }

    if (any(y <= 0) || any(sc <= 0))
      return(10^6)

    if (abs(xi) >= 10^(-5)) {
      ooxi <- 1 / xi
    } else {
      ooxi <- sign(xi) * 10^5
    }

    zz <- nsample * (log(sc)) + sum(exp(ooxi * log(y))) + sum(log(y) * (1 - (ooxi)))

    # Coles-Dixon penalty
    zz <- zz - log(cd.hos(xi))

    return(zz)
  }

  mup <- abs(mean(xdat)) * 3
  sigup <- sqrt(var(xdat)) * 5

  tryCatch(
    for (i in 1:nrow(init)) {
      value <- try(solnp(init[i, ], fun = gev1.lik.CD.h,
                         LB = c(-mup, 0, -0.99), UB = c(mup, sigup, 0.499),
                         control = list(trace = 0, outer.iter = 20,
                                        delta = 1.e-7, inner.iter = 200, tol = 1.e-5)))

      if (is(value)[1] == "try-error") {
        k[[i]] <- list(value = 10^6)
      } else {
        k[[i]] <- value
      }
    }
  )

  optim_value <- data.frame(num = 1:n, value = sapply(k, function(x) x$value[which.min(x$value)]))
  optim_table1 <- optim_value[order(optim_value$value), ]
  selc_num <- optim_table1[1, "num"]

  x <- k[[selc_num]]

  z$conv <- x$convergence
  z$nllh <- x$value[which.min(x$value)]
  z$mle <- x$par  # Hosking style parameter

  return(z)
}


#' Coles-Dixon Penalty Function
#'
#' Computes the Coles-Dixon prior probability for the shape parameter.
#'
#' @param sxi Shape parameter value (Hosking style, negative for heavy tails).
#'
#' @return Prior probability value.
#'
#' @keywords internal
cd.hos <- function(sxi) {
  if (sxi >= 0) {
    pf <- 1.0
  } else if (sxi <= -1) {
    pf <- 1.e-20
  } else {
    pf <- exp(-((1 / (1 + sxi)) - 1))  # Hosking style xi
  }
  return(pf)
}


#' Restricted MLE for GEV (Mixed Estimation)
#'
#' Computes restricted maximum likelihood estimates for GEV parameters
#' with constraints on the mean or median matching the sample statistics.
#'
#' @param xdat Numeric vector of data.
#' @param ntry Number of random starting points. Default is 5.
#' @param rest Restriction type: 'mean' (default) or 'median'.
#' @param quant Probabilities for quantile estimation. Default is c(0.99, 0.995).
#' @param trim Left trimming level. Default is 0.
#' @param CD.mle Coles-Dixon MLE (optional). If NULL, computed internally.
#' @param mle Standard MLE (optional). If NULL, computed internally.
#' @param second Logical. If TRUE, compute second-stage REMLE. Default is TRUE.
#' @param w.mpse Logical. If TRUE, compute MPSE. Default is FALSE.
#'
#' @return A list containing:
#'   \item{remle1}{First-stage REMLE estimates}
#'   \item{qua.remle1}{Quantiles from first-stage REMLE}
#'   \item{remle2}{Second-stage REMLE estimates (if second=TRUE)}
#'   \item{qua.remle2}{Quantiles from second-stage REMLE}
#'   \item{rest.method}{Restriction method used}
#'
#' @keywords internal
gev.remle <- function(xdat, ntry = 5, rest = 'mean', quant = c(0.99, 0.995),
                      trim = 0, CD.mle = NULL, mle = NULL, second = TRUE,
                      w.mpse = FALSE) {
  z <- list()
  k <- list()
  data <- xdat
  if (rest == 'med') rest <- 'median'

  nsample <- length(xdat)
  z$nsample <- nsample

  if (is.null(CD.mle)) {
    CD.mle <- rep(NA, 3)
    CD.mle <- gev1.CD(xdat = data, ntry = 5)$mle  # Hosking style xi
  }
  if (is.null(mle)) {
    mle <- rep(NA, 3)
    mle <- gev.max(xdat = data, ntry = 5)$mle  # Hosking style xi
  }

  if (ntry < 4) ntry <- 4
  init <- matrix(NA, ntry, 3)
  init <- ginit.max(data = xdat, ntry = ntry - 2)
  init <- rbind(init, matrix(NA, 2, 3))
  init[ntry - 1, 1:3] <- mle[1:3]
  init[ntry, 1:3] <- CD.mle[1:3]

  init <- na.omit(init)
  ntry <- nrow(init)

  upsig <- sqrt(var(xdat)) * 3
  upmu <- abs(mean(xdat)) * 3

  ntrim <- trim
  med.dat <- median(data)
  mean.dat <- mean(data)
  lmx <- lmoms(data, nmom = 3)

  nar <- c('rest.Ex', 'rest.med')
  if (rest == 'mean') {
    sel <- 1
  } else if (rest == 'median') {
    sel <- 2
  }

  # Define constraint functions in parent environment
  rest.Ex <- function(a, xdat, lmx, med.dat, mean.dat) {
    lam1 <- lmomgev(para = vec2par(a, 'gev'))$lambdas[1]
    return(lam1 - mean.dat)
  }

  rest.med <- function(a, xdat, lmx, med.dat, mean.dat) {
    med <- quagev(0.5, vec2par(a, 'gev'))
    return(med - med.dat)
  }

  rest.Ex2 <- function(a, xdat, lmx, med.dat, mean.dat) {
    lamp <- lmomgev(para = vec2par(a, 'gev'))$lambdas[1:2]
    return(c(lamp[1] - lmx$lambdas[1], lamp[2] - lmx$lambdas[2]))
  }

  rest.med2 <- function(a, xdat, lmx, med.dat, mean.dat) {
    medx <- quagev(0.5, vec2par(a, 'gev'))
    lam2 <- lmomgev(para = vec2par(a, 'gev'))$lambdas[2]
    return(c(medx - med.dat, lam2 - lmx$lambdas[2]))
  }

  gev.lik.remax <- function(a, xdat, lmx, med.dat, mean.dat) {
    mu <- a[1]
    sc <- a[2]
    xi <- a[3]

    y <- (xdat - mu) / sc
    y <- 1 - xi * y  # Hosking style

    nsample <- length(xdat)
    for (i in 1:nsample) {
      y[i] <- max(0, y[i], na.rm = TRUE)
    }

    if (any(y <= 0) || any(sc <= 0))
      return(10^6)

    if (abs(xi) >= 10^(-5)) {
      ooxi <- 1 / xi
    } else {
      ooxi <- sign(xi) * 10^5
    }

    zz <- nsample * (log(sc)) + sum(exp(ooxi * log(y))) + sum(log(y) * (1 - (ooxi)))
    return(zz)
  }

  # Select constraint function
  if (sel == 1) {
    eq_fun <- rest.Ex
  } else {
    eq_fun <- rest.med
  }

  tryCatch(
    for (i in 1:ntry) {
      work <- list()
      work <- try(solnp(init[i, ], fun = gev.lik.remax,
                        eqfun = eq_fun, eqB = 0,
                        LB = c(-upmu, 0, -0.99), UB = c(upmu, upsig, 0.99),
                        control = list(trace = 0, outer.iter = 20,
                                        delta = 1.e-6, inner.iter = 200, tol = 1.e-5),
                        xdat = xdat, lmx = lmx, med.dat = med.dat,
                        mean.dat = mean.dat))

      if (is(work)[1] == "try-error") {
        k[[i]] <- list(values = 10^6)
      } else if (work$convergence != 0) {
        k[[i]] <- list(values = 10^6)
      } else {
        k[[i]] <- work
      }
    }
  )

  optim_value <- data.frame(num = 1:ntry, values = sapply(k, function(x) x$values[which.min(x$values)]))
  optim_table1 <- optim_value[order(optim_value$values), ]
  selc_num <- optim_table1[1, "num"]
  x <- list()

  x <- k[[selc_num]]

  z$remle1.value <- x$values[which.min(x$values)]
  z$remle1 <- x$pars

  if (is.null(x$pars) | any(is.na(x$pars)) | z$remle1.value == 10^6) {
    x$pars <- NA
    z$remle1 <- NA
    z$qua.remle1 <- NA
  } else if (x$pars[2] <= 0 | abs(x$pars[3]) >= 1.0) {
    x$pars <- NA
    z$remle1 <- NA
    z$qua.remle1 <- NA
  } else {
    z$qua.remle1 <- quagev(quant, vec2par(x$pars, 'gev'))
  }

  z$rest.method <- rest

  if (second == TRUE) {
    k2 <- list()
    if (sel == 1) {
      eq_fun2 <- rest.Ex2
    } else {
      eq_fun2 <- rest.med2
    }

    init2 <- rbind(init, rep(NA, 3))
    init2[nrow(init2), 1:3] <- z$remle1[1:3]
    init2 <- na.omit(init2)
    xtry <- nrow(init2)

    tryCatch(
      for (i in 1:xtry) {
        value2 <- list()
        value2 <- try(solnp(init2[i, ], fun = gev.lik.remax,
                            eqfun = eq_fun2, eqB = c(0, 0),
                            LB = c(-upmu, 0, -0.99), UB = c(upmu, upsig, 0.99),
                            control = list(trace = 0, outer.iter = 20,
                                            delta = 1.e-6, inner.iter = 200, tol = 1.e-5),
                            xdat = xdat, lmx = lmx, med.dat = med.dat,
                            mean.dat = mean.dat))

        if (is(value2)[1] == "try-error") {
          k2[[i]] <- list(values = 10^6)
        } else if (value2$convergence != 0) {
          k2[[i]] <- list(values = 10^6)
        } else {
          k2[[i]] <- value2
        }
      }
    )

    optim_value <- data.frame(num = 1:xtry, values = sapply(k2, function(x) x$values[which.min(x$values)]))
    optim_table1 <- optim_value[order(optim_value$values), ]
    selc_num <- optim_table1[1, "num"]
    xw <- list()

    xw <- k2[[selc_num]]

    z$remle2.value <- xw$values[which.min(xw$values)]
    z$remle2 <- xw$pars

    if (is.null(xw$pars) | any(is.na(xw$pars)) | z$remle2.value == 10^6) {
      xw$pars <- NA
      z$remle2 <- NA
      z$qua.remle2 <- NA
    } else if (xw$pars[2] <= 0 | abs(xw$pars[3]) >= 1.0) {
      xw$pars <- NA
      z$remle2 <- NA
      z$qua.remle2 <- NA
    } else {
      z$qua.remle2 <- quagev(quant, vec2par(xw$pars, 'gev'))
    }
  }

  z$mpse <- NA
  z$qua.mpse <- NA

  return(z)
}


# ============================================================
# Weight computation functions (from weight.com.BMA_15Dec25.R)
# ============================================================

#' @keywords internal
prior.beta.stnary.ma <- function(x = NULL, xi_lme = NULL,
                                  p = 2, q = 5, al, bl) {
  Bef <- function(w) { ((-al + w)^(p - 1)) * ((bl - w)^(q - 1)) }
  Be <- integrate(Bef, lower = al, upper = bl)[1]$value

  if ((-al + x) <= 0) {
    left.w <- 0
  } else {
    left.w <- (-al + x)
  }
  if ((bl - x) <= 0) {
    right.w <- 0
  } else {
    right.w <- (bl - x)
  }

  prior <- (left.w^(p - 1)) * (right.w^(q - 1)) / Be
  return(prior)
}


#' @keywords internal
set.prior <- function(pen = NULL, numk = NULL, xi_lme = NULL, kpar = NULL,
                      weight = NULL) {
  zp <- list()
  prior <- rep(NA, numk)
  lme <- xi_lme

  if (weight == 'like') {
    if (pen == "beta") {
      p <- 2
      c0 <- 0.35
      c1 <- 40
      c2 <- 12
    } else if (pen == "norm") {
      multa <- 2.1
      lme.cut <- -0.5
      std.cut <- -0.45
    }
  } else if (weight == 'gLd' | weight == 'med') {
    if (pen == "beta") {
      p <- 2
      c0 <- 0.35
      c1 <- 20
      c2 <- 7
    } else if (pen == "norm") {
      multa <- 1.5
      lme.cut <- -0.45
      std.cut <- -0.4
    }
  }

  if (pen == "norm") {
    if (lme >= 0) {
      prior <- rep(1 / numk, numk)
      mu.prior <- 0
      std.prior <- 1
    }

    if (lme < 0) {
      if (lme >= lme.cut) {
        mu.prior <- multa * lme
      } else if (lme < lme.cut) {
        mu.prior <- multa * lme.cut
      }

      if (weight == "like") {
        if (lme >= std.cut) {
          std.prior <- -(std.cut - lme) / 5 + 0.11
        } else {
          std.prior <- 0.11
        }
      }
      if (weight == "gLd" | weight == "med") {
        if (lme >= std.cut) {
          std.prior <- -(std.cut - lme) / 4 + 0.14
        } else {
          std.prior <- 0.14
        }
      }

      for (i in 1:numk) {
        prior[i] <- dnorm(x = kpar[i], mean = mu.prior, sd = std.prior)
      }
    }
  } else if (pen == "beta") {
    aa <- max(-1, xi_lme - c0)
    bb <- min(0.5, xi_lme + c0)
    al <- min(aa, bb)
    bl <- max(aa, bb)

    if (xi_lme <= 0) {
      qlim <- min(abs(xi_lme) * c1, c2)
    } else {
      qlim <- 0.0
    }

    p <- p
    q <- p + qlim

    for (i in 1:numk) {
      prior[i] <- prior.beta.stnary.ma(x = kpar[i], xi_lme = lme,
                                        p = p, q = q, al = al, bl = bl)
    }
  }

  zp$prior <- prior
  if (pen == "norm") zp$prior_mu_std <- c(mu.prior, std.prior)
  if (pen == "beta") zp$p_q_beta <- round(c(p, q), 2)
  return(zp)
}


#' @keywords internal
weight.com.new <- function(data = NULL, numk = NULL, hosking = NULL,
                           kpar = NULL, numom = 3, xqa = NULL, varcom = TRUE,
                           boot.lme = TRUE, cov.lme = NULL, surr = FALSE,
                           type = 'full', trim = NULL, cov.type = 'ratio',
                           bma = TRUE, pen = "norm") {
  z <- list()
  quant <- hosking$quant
  weight <- hosking$weight
  B <- hosking$B
  start <- hosking$start

  if (bma == TRUE) {
    prior <- set.prior(pen = pen, numk = numk, xi_lme = hosking$lme[3],
                       kpar = kpar, weight = weight)$prior
  } else if (bma != TRUE) {
    prior <- rep(1 / numk, numk)
  }

  prob.call <- list()
  prob.call <- dist.noboot(data = data, numk = numk, hosking = hosking,
                           boot.lme = boot.lme,
                           kpar = kpar, numom = numom, ntry = 5, varcom = varcom,
                           cov.lme = cov.lme, trim = trim, cov.type = cov.type)

  notid <- which(!is.na(prob.call$mle3[, 1]))
  prob.call$mle3 <- prob.call$mle3[notid, ]
  mle3 <- prob.call$mle3
  numk <- length(notid)

  kpar <- kpar[notid]

  if (surr == TRUE) {
    nxq <- length(xqa)
    xzp <- matrix(NA, nrow = nxq, ncol = numk)
  }
  numq <- length(quant)

  zp <- matrix(NA, nrow = numq, ncol = numk)
  wtgd <- rep(NA, numk)
  bmaw <- rep(NA, numk)

  lme <- hosking$lme
  mle <- hosking$mle

  if (weight == 'gLd' | weight == 'med' | weight == 'cvt' | weight == 'fcv') {
    if (varcom == TRUE | weight == 'fcv') {
      cov2 <- list()
      for (i in 1:numk) cov2[[i]] <- prob.call$cov2[[notid[i]]]
    }

    if (weight == 'cvt') {
      naic <- prob.call$aic[notid]
      id <- which(is.na(naic))
      amin <- naic - min(naic, na.rm = TRUE)
      wtgd <- exp(-amin) / sum(exp(-amin), na.rm = TRUE)
      wtgd[id] <- 0
      bmaw <- prior * exp(-amin) / sum(prior * exp(-amin), na.rm = TRUE)
      bmaw[id] <- 0
    }

    if (weight == 'gLd' | weight == 'med' | weight == 'fcv') {
      if (weight == 'gLd') {
        kcol <- 1
      } else if (weight == 'med') {
        kcol <- 2
      }
      if (weight == 'fcv') {
        kcol <- 3
      }

      prob.mtx <- prob.call$prob.mtx[notid, ]
      wtgd[1:numk] <- prob.mtx[1:numk, kcol] / sum(prob.mtx[, kcol], na.rm = TRUE)
      bmaw[1:numk] <- prior * prob.mtx[1:numk, kcol] / sum(prior * prob.mtx[, kcol], na.rm = TRUE)
    }
  }

  if (weight == 'lcv' | weight == 'like' | weight == 'opt') {
    pertr <- hosking$pertr

    if (weight == 'lcv' | weight == 'like') {
      prob.call <- wlik.xifix(data, numk = numk, kpar = kpar, weight = weight,
                              pertr = pertr, varcom = varcom, type = type,
                              prior = prior, trim = trim)

      if (varcom == TRUE) {
        cov2 <- list()
        for (i in 1:numk) cov2[[i]] <- prob.call$cov2[[notid[i]]]
      }
    }

    wtgd <- prob.call$wt
    bmaw <- prob.call$bmaw
    prob.call$mle3 <- t(prob.call$kfix)

    for (ip in 1:numk) {
      save <- vec2par(prob.call$kfix[1:3, ip], 'gev')
      zp[1:numq, ip] <- quagev(quant[1:numq], save)
      if (surr == TRUE) xzp[1:nxq, ip] <- quagev(xqa[1:nxq], save)
    }
  } else {
    for (ip in 1:numk) {
      if (any(is.na(prob.call$mle3[ip, 1:3]))) {
        zp[1:numq, ip] <- NA
        if (surr == TRUE) xzp[1:nxq, ip] <- NA
      } else {
        save <- vec2par(prob.call$mle3[ip, 1:3], 'gev')
        zp[1:numq, ip] <- quagev(quant[1:numq], save)
        if (surr == TRUE) xzp[1:nxq, ip] <- quagev(xqa[1:nxq], save)
      }
    }
  }

  z$weight <- weight
  z$numk <- numk
  z$wtgd <- wtgd
  z$bmaw <- bmaw
  z$zp <- zp
  z$kpar <- kpar
  z$prob.call <- prob.call

  if (varcom == TRUE) {
    z$prob.call$cov2 <- cov2
  } else {
    z$prob.call$cov2 <- NA
  }
  if (surr == TRUE) z$xzp <- xzp

  return(z)
}


#' @keywords internal
wlik.xifix <- function(data = NULL, numk = NULL, kpar = NULL, weight = NULL,
                       pertr = 0.85, varcom = FALSE, type = NULL,
                       prior = NULL, trim = NULL) {
  z <- list()
  eprior <- prior
  kfix <- matrix(NA, nrow = 3, ncol = numk)
  nbic <- rep(NA, numk)
  bmin <- rep(NA, numk)
  wt <- rep(NA, numk)
  bmaw <- rep(NA, numk)
  para <- rep(NA, 3)
  lm3 <- list()
  re <- sort(data)
  nsam <- length(data)

  ntrim <- trim
  lmtrim <- lmoms(re[(ntrim + 1):nsam], nmom = 3)
  if (is.null(eprior)) eprior <- 1

  cov2.cv <- list()

  for (ip in 1:numk) {
    kfix[1:3, ip] <- pargev.xifix.ma(lmtrim, xifix = kpar[ip])$para[1:3]

    if (varcom == TRUE) {
      Hess <- PrescottW(par2 = as.vector(kfix[1:2, ip]), xifix = kpar[ip],
                        nsam = length(re[(ntrim + 1):nsam]))

      if (is.na(Hess[1, 1])) {
        cov2.cv[[ip]] <- matrix(NA, 2, 2)
      } else if (det(Hess) <= 0) {
        cov2.cv[[ip]] <- matrix(NA, 2, 2)
      } else {
        cov2.cv[[ip]] <- solve(Hess)
      }
    }
  }

  if (weight == 'like') {
    for (ip in 1:numk) {
      ntrim <- trim
      dtrim <- re[(ntrim + 1):nsam]
      nbic[ip] <- gev.xilik(a = kfix[1:2, ip], xifix = kpar[ip], xdat = dtrim)
      if (nbic[ip] >= 10^6) nbic[ip] <- NA
    }

    numpar <- 2
    nbic <- 2 * nbic + 2 * numpar
    bmin <- nbic - min(nbic, na.rm = TRUE)
    wt <- exp(-bmin) / sum(exp(-bmin), na.rm = TRUE)
    bmaw <- prior * exp(-bmin) / sum(prior * exp(-bmin), na.rm = TRUE)
  }

  bmaw[which(is.na(bmaw))] <- 0
  wt[which(is.na(wt))] <- 0

  z$bmaw <- bmaw
  z$wt <- wt
  z$aic <- nbic
  z$kfix <- kfix
  if (varcom == TRUE) {
    z$cov2 <- cov2.cv
  }

  return(z)
}


#' @keywords internal
dist.noboot <- function(data = NULL, numk = NULL, hosking = NULL, boot.lme = TRUE,
                        kpar = NULL, numom = NULL, ntry = 5, varcom = NULL,
                        cov.lme = NULL, trim = NULL, cov.type = 'lambda') {
  cov2 <- list()
  kfix <- list()
  zw <- list()
  weight <- hosking$weight
  B <- hosking$B
  if (boot.lme == TRUE) boot.lme <- TRUE

  mle3 <- matrix(NA, nrow = numk, ncol = 3)
  lmB.med <- rep(NA, B)
  aic <- rep(NA, numk)

  nsample <- length(data)
  lm <- lmoms(data, nmom = numom)
  medata <- median(data)

  for (ip in 1:numk) {
    if (is.na(kpar[ip])) {
      aic[ip] <- NA
      mle3[ip, 1:3] <- NA
    } else {
      kfix[[ip]] <- gev.xifix.sing(xdat = data, xifix = kpar[ip], ntry = ntry,
                                   varcom = varcom)

      if (kfix[[ip]]$conv != 0) {
        aic[ip] <- NA
        mle3[ip, 1:3] <- NA
      } else {
        aic[ip] <- 2 * kfix[[ip]]$nllh + 2 * 2
        mle3[ip, 1:3] <- kfix[[ip]]$mle[1:3]
      }
      if (varcom == TRUE) cov2[[ip]] <- kfix[[ip]]$cov
    }
  }

  if (all(is.na(mle3[1:numk, 1]))) {
    lm3 <- lmoms(data, nmom = 3)
    for (ip in 1:numk) {
      mle3[ip, 1:3] <- pargev.xifix.ma(lm3, xifix = kpar[ip])$para[1:3]
    }
  }

  if (weight == 'gLd') {
    if (cov.type == 'ratio') {
      if (boot.lme == TRUE) {
        cov.lm <- hosking$cov.rat
      } else {
        cov.lm <- cov.lme$cov.rat
      }
    } else if (cov.type == 'lambda') {
      if (boot.lme == TRUE) {
        cov.lm <- hosking$cov.lambda
      } else {
        cov.lm <- cov.lme$cov.lambda
      }
    }
    Vinv <- solve(cov.lm)
    detV <- det(cov.lm)
    zw$cov.lm <- list(Vinv = Vinv, detV = detV)
  }

  if (weight == 'med') {
    if (cov.type == 'ratio') {
      if (boot.lme == TRUE) {
        cov.med <- hosking$cov.med.rat
      } else {
        cov.med <- cov.lme$cov.med.rat
      }
    } else if (cov.type == 'lambda') {
      if (boot.lme == TRUE) {
        cov.med <- hosking$cov.med.lam
      } else {
        cov.med <- cov.lme$cov.med.lam
      }
    }

    Vinv.med <- solve(cov.med)
    detV.med <- det(cov.med)
    Vinv <- Vinv.med
    detV <- detV.med
    zw$cov.lm <- list(Vinv = Vinv.med, detV = detV.med)
  }

  if (weight == 'gLd' | weight == 'med') {
    normal <- list()
    normal$prob.mtx <- matrix(NA, nrow = numk, ncol = 2)

    normal <- com.prdist(data = data, numk = numk, kfix = kfix,
                         Vinv = Vinv, detV = detV, numom = numom,
                         hosking = hosking, trim = trim,
                         cov.type = cov.type)

    for (ip in 1:numk) {
      if (any(is.na(kfix[[ip]]$mle[1:3]))) {
        mle3[ip, 1:3] <- NA
        normal$prob.mtx[ip, 1:2] <- 0.0
      } else {
        mle3[ip, 1:3] <- kfix[[ip]]$mle[1:3]
      }
    }
  }

  zw$aic <- aic
  zw$mle3 <- mle3
  zw$kfix <- kfix
  if (varcom == TRUE) {
    zw$cov2 <- cov2
  } else {
    zw$cov2 <- NA
  }
  zw$prob.mtx <- matrix(NA, numk, 3)

  if (weight == 'gLd' | weight == 'med') {
    zw$prob.mtx <- normal$prob.mtx
    zw$gdd <- normal$gdd
  }

  return(zw)
}


#' @keywords internal
com.prdist <- function(data = NULL, numk = NULL, kfix = NULL,
                       Vinv = NULL, detV = NULL, numom = NULL,
                       hosking = NULL, trim = NULL, cov.type = "lambda") {
  weight <- hosking$weight

  ww1 <- list()
  kfix.lm <- list()
  dd <- matrix(NA, nrow = numom, ncol = 2)
  gdd <- matrix(NA, nrow = numk, ncol = 2)
  simple <- rep(NA, numk)
  prob.mtx <- matrix(NA, nrow = numk, ncol = 2)

  re <- sort(data)
  nsam <- length(data)
  ntrim <- trim

  lm <- lmoms(data, nmom = numom)

  for (ip in 1:numk) {
    if (kfix[[ip]]$conv != 0) {
      prob.mtx[ip, 1:2] <- 0
    } else {
      if (kfix[[ip]]$mle[3] <= -1.0) kfix[[ip]]$mle[3] <- -0.999
      simple[ip] <- kfix[[ip]]$mle[3]
      savek <- vec2par(kfix[[ip]]$mle, 'gev')
      kfix.lm[[ip]] <- lmomgev(savek)
    }
  }

  if (weight == 'gLd') {
    for (ip in 1:numk) {
      if (kfix[[ip]]$conv == 0) {
        lmtrim.cvt <- lmoms(re[(ntrim + 1):nsam], nmom = 3)
        dd[1, 1] <- lmtrim.cvt$lambdas[1] - kfix.lm[[ip]]$lambdas[1]

        if (cov.type == 'ratio') {
          dd[2:numom, 1] <- lmtrim.cvt$ratios[2:numom] - kfix.lm[[ip]]$ratios[2:numom]
        } else if (cov.type == 'lambda') {
          dd[2, 1] <- lmtrim.cvt$lambdas[2] - kfix.lm[[ip]]$lambdas[2]
          dd[3, 1] <- lmtrim.cvt$lambdas[3] - kfix.lm[[ip]]$lambdas[3]
        }

        gdd[ip, 1] <- t(dd[1:numom, 1]) %*% Vinv %*% dd[1:numom, 1]
        prob.mtx[ip, 1] <- exp(-gdd[ip, 1] / 2) / ((2 * pi)^(numom / 2) * sqrt(detV))
        if (is.na(prob.mtx[ip, 1])) prob.mtx[ip, 1] <- 0
      }
    }
  }

  if (weight == 'med') {
    for (ip in 1:numk) {
      if (kfix[[ip]]$conv == 0) {
        medata <- median(re[(ntrim + 1):nsam])
        lmtrim.cvt <- lmoms(re[(ntrim + 1):nsam], nmom = 3)

        savek <- vec2par(kfix[[ip]]$mle, 'gev')
        dd[1, 2] <- medata - quagev(0.5, savek)

        if (cov.type == 'ratio') {
          dd[2:numom, 2] <- lmtrim.cvt$ratios[2:numom] - kfix.lm[[ip]]$ratios[2:numom]
        } else if (cov.type == 'lambda') {
          dd[2, 2] <- lmtrim.cvt$lambdas[2] - kfix.lm[[ip]]$lambdas[2]
          dd[3, 2] <- lmtrim.cvt$lambdas[3] - kfix.lm[[ip]]$lambdas[3]
        }

        gdd[ip, 2] <- t(dd[1:numom, 2]) %*% Vinv %*% dd[1:numom, 2]
        prob.mtx[ip, 2] <- exp(-gdd[ip, 2] / 2) / ((2 * pi)^(numom / 2) * sqrt(detV))
        if (is.na(prob.mtx[ip, 2])) prob.mtx[ip, 2] <- 0
      }
    }
  }

  ww1$prob.mtx <- prob.mtx
  ww1$gdd <- gdd
  return(ww1)
}


#' @keywords internal
new.kpar2 <- function(wtgd = NULL, numk = NULL, kpar = NULL,
                      remove = 0.004, dist = 0.02, tre1 = 0.01, tre8 = 0.5) {
  wtgd[which(is.na(wtgd))] <- 0
  zp <- list()
  aw <- 0

  nid <- which(wtgd < remove)
  id <- which(wtgd >= remove)
  kpar <- kpar[id]
  wtgd <- wtgd[id]
  numk2 <- length(kpar)
  if (numk2 < numk) aw <- 10
  numk <- numk2
  kpar2 <- kpar

  if (numk == 0) {
    kpar2 <- -0.1
    numk <- 1
  }

  if (numk == 1) {
    ad <- 1.5
    kpar2 <- c(NA, NA, kpar2, NA, NA)
    kpar2[2] <- max(kpar2[3] - dist / ad, -0.99)
    kpar2[1] <- max(kpar2[3] - 2 * dist / ad, -0.99)
    kpar2[4] <- min(kpar2[3] + dist / ad, 0.49)
    kpar2[5] <- min(kpar2[3] + 2 * dist / ad, 0.49)
    aw <- 11
  } else {
    work1 <- work2 <- work3 <- 0

    if (wtgd[1] >= tre1 & wtgd[1] <= tre8) {
      kpar2 <- c(rep(NA, 2), kpar)
      kpar2[2] <- max(kpar[1] - dist, -0.99)
      kpar2[1] <- max(kpar2[2] - dist, -0.99)
      work1 <- 1
    } else if (wtgd[1] > tre8) {
      kpar2 <- c(rep(NA, 4), kpar)
      kpar2[4] <- max(kpar[1] - dist / 2, -0.99)
      kpar2[3] <- max(kpar2[4] - dist / 2, -0.99)
      kpar2[2] <- max(kpar2[3] - dist / 2, -0.99)
      kpar2[1] <- max(kpar2[2] - dist / 2, -0.99)
      work1 <- 20
    }

    if (work1 == 0 & wtgd[numk] >= tre1 & wtgd[numk] <= tre8) {
      kpar2 <- c(kpar, rep(NA, 2))
      kpar2[numk + 1] <- min(kpar2[numk] + dist, 0.49)
      kpar2[numk + 2] <- min(kpar2[numk + 1] + dist, 0.49)
      work2 <- 1
    }
    if (work1 == 0 & wtgd[numk] > tre8) {
      kpar2 <- c(kpar, rep(NA, 4))
      kpar2[numk + 1] <- min(kpar2[numk] + dist / 2, 0.49)
      kpar2[numk + 2] <- min(kpar2[numk + 1] + dist / 2, 0.49)
      kpar2[numk + 3] <- min(kpar2[numk + 2] + dist / 2, 0.49)
      kpar2[numk + 4] <- min(kpar2[numk + 3] + dist / 2, 0.49)
      work2 <- 20
    }

    if (work1 == 1 & wtgd[numk] >= tre1) {
      kpar2 <- c(kpar2, rep(NA, 2))
      kpar2[numk + 3] <- min(kpar2[numk + 2] + dist, 0.49)
      kpar2[numk + 4] <- min(kpar2[numk + 3] + dist, 0.49)
      work2 <- 2
    }

    lek <- length(kpar2)
    mxid <- which.max(wtgd)
    if (wtgd[mxid] > tre8 & mxid != 1 & mxid != numk) {
      kpar2 <- c(kpar2, rep(NA, 4))
      kpar2[lek + 1] <- max(kpar[mxid] - dist / 3, -0.99)
      kpar2[lek + 2] <- max(kpar[mxid] - 2 * dist / 3, -0.99)
      kpar2[lek + 3] <- min(kpar[mxid] + dist / 3, 0.49)
      kpar2[lek + 4] <- min(kpar[mxid] + 2 * dist / 3, 0.49)
      work3 <- 4
    }
    aw <- aw + work1 + work2 + work3
  }

  numk <- length(kpar2)

  zp$kpar2 <- sort(kpar2)
  zp$aw <- aw
  zp$numk <- numk
  return(zp)
}


# ============================================================
# Candidate xi selection functions (from cand.xi.paper_12Dec25.R)
# ============================================================

#' @keywords internal
cand.xi.new.paper <- function(data, hosking = NULL, mle = NULL, pick0 = 0.95,
                              nint = 256, start = 'mle', numk = NULL,
                              figure = TRUE, cov.lme = NULL, bma = FALSE,
                              pen = "beta") {
  wx <- list()
  SMALL <- 1e-05
  kpar <- rep(NA, numk)
  move <- FALSE

  while (start == 'mle' | start == 'mix') {
    newse <- max(0.3, abs(mle[3])) * 0.2

    mle.coles <- mle
    mle.coles[3] <- -mle[3]

    xlow <- mle.coles[3] - 25 * newse
    xup <- mle.coles[3] + 12 * newse
    xlow <- max(-0.99, min(0.5, xlow), na.rm = TRUE)
    xup <- min(0.99, max(-0.5, xup), na.rm = TRUE)
    xlow2 <- min(xlow, xup)
    xup2 <- max(xlow, xup)

    pilo <- (1 - pick0) / 2
    piup <- 1 - pilo
    eqpr <- seq(pilo, piup, by = pick0 / (numk - 1))
    numpi <- floor(numk / 2)
    pipick <- rep(NA, numpi)

    pipick[1] <- pick0
    for (i in 2:numpi) {
      pipick[i] <- 1 - 2 * eqpr[i]
    }

    get.ci <- gev.profxi.mdfy.paper(data, mle = mle.coles,
                                    xlow = xlow2, xup = xup2,
                                    pick.v = pipick, nint = nint, figure = figure)

    if (get.ci$fail == TRUE) {
      start <- 'lme'
      move <- TRUE
      break
    }

    wx$get.ci <- get.ci

    lo.lim <- min(-get.ci$w1$ci1[1], -get.ci$w1$ci2[1])
    up.lim <- max(-get.ci$w1$ci1[1], -get.ci$w1$ci2[1])

    kpar[1] <- max(lo.lim, -0.99)
    kpar[numk] <- min(up.lim, 0.99)

    if (is.odd.me(numk)) kpar[floor(numk / 2) + 1] <- -get.ci$w1$xmax

    for (i in 2:numpi) {
      lo.lim <- min(-get.ci$w1$ci1[i], -get.ci$w1$ci2[i])
      up.lim <- max(-get.ci$w1$ci1[i], -get.ci$w1$ci2[i])
      kpar[i] <- max(lo.lim, -0.99)
      kpar[numk - i + 1] <- min(up.lim, 0.99)
    }

    if (kpar[numk] < -0.3) {
      kpar2 <- rep(NA, numk + 2)
      kpar2[numk + 1] <- max(-0.3, kpar[numk] + 0.05)
      kpar2[numk + 2] <- max(-0.3, kpar2[numk + 1] + 0.05)
      kpar[1:numk] <- c(kpar[3:numk], kpar2[numk + 1], kpar2[numk + 2])
    }

    wx$kpar.mle <- kpar
    break
  }

  if (start == 'lme' | start == 'mix') {
    pilo <- (1 - pick0) / 2
    piup <- 1 - pilo
    eqpr <- seq(pilo, piup, by = pick0 / (numk - 1))

    if (!is.null(cov.lme)) {
      hosking$lme.boot <- cov.lme$lme.boot
    } else {
      if (hosking$boot.lme == FALSE & move == TRUE) {
        Bnew <- hosking$B
        qua <- hosking$quant
        hosking <- lme.boots.new(data = data, B = Bnew, quant = qua, boot = TRUE)
      }
    }

    kpar <- quantile(hosking$lme.boot[, 3], prob = eqpr, na.rm = TRUE)

    if (kpar[numk] < -0.5) kpar[numk] <- max(-0.5, kpar[numk] + 0.05)
    wx$kpar.lme <- kpar
    if (move == TRUE) wx$kpar.mle <- kpar
    hosking$start <- start
  }

  if (start == 'mix') {
    wx$kpar.mix <- (wx$kpar.mle + wx$kpar.lme) / 2
  }

  if (start == 'mle') {
    kpar <- wx$kpar.mle
  } else if (start == 'lme') {
    kpar <- wx$kpar.lme
  } else if (start == 'mix') {
    kpar <- wx$kpar.mix
  }

  if (any(abs(kpar) <= SMALL)) {
    kpar[which(abs(kpar) <= SMALL)] <- SMALL * 1.2
  }

  if (kpar[numk] >= 0.5) kpar[numk] <- 0.49
  wx$kpar <- kpar
  wx$start <- start
  if (start == 'mle') wx$ymin <- get.ci$ymin

  return(wx)
}


#' @keywords internal
is.odd.me <- function(x) {
  return(ifelse(x %% 2 == 1, TRUE, FALSE))
}


#' @keywords internal
gev.profxi.mdfy.paper <- function(data = NULL, mle = NULL, xlow, xup,
                                  pick.v = NULL, nint = 256,
                                  figure = FALSE) {
  pick <- pick.v[1]
  w <- list()
  v <- numeric(nint)

  w$fail <- FALSE
  w$start <- 'mle'

  conf <- pick.v
  xlow.t <- min(xlow, xup)
  xup.t <- max(xlow, xup)
  xlow <- max(-0.9, xlow.t, na.rm = TRUE)
  xup <- min(0.9, xup.t, na.rm = TRUE)
  x <- seq(xlow, xup, length.out = nint)
  sol <- c(mle[1], mle[2])

  gev.plikxi <- function(a, data = data, xi = xi) {
    if (a[2] <= 0) {
      l <- 10^6
      return(l)
    }
    if (abs(xi) < 10^(-6)) {
      y <- (data - a[1]) / a[2]
      l <- length(y) * log(a[2]) + sum(exp(-y)) + sum(y)
    } else {
      y <- (data - a[1]) / a[2]
      y <- 1 + xi * y
      if (any(y <= 0)) {
        l <- 10^6
      } else {
        l <- length(y) * log(a[2]) + sum(y^(-1 / xi)) + sum(log(y)) * (1 / xi + 1)
      }
      l
    }
  }

  for (i in 1:nint) {
    xi <- x[i]
    opt <- optim(sol, gev.plikxi, data = data, xi = xi)
    sol <- opt$par
    v[i] <- opt$value
  }

  for (i in 1:nint) {
    if (v[i] >= 10^6) v[i] <- NA
  }

  d <- data.frame(x = x, v = -v)
  d <- na.omit(d)

  irt <- 0
  idmax <- which.max(d$v)
  if (idmax == 1) {
    irt <- 1
  } else if (idmax == length(d$x)) {
    irt <- 1
  }

  if (irt == 1) {
    w$start <- 'lme'
    w$fail <- TRUE
    return(w)
  }

  port <- 1.3
  x0 <- min(d$x)
  x9 <- max(d$x)
  bd0 <- max(0.2, abs(x0)) * sign(x0)
  bd9 <- max(0.2, abs(x9)) * sign(x9)

  xmin <- min(-0.7, x0 - sign(x0) * port * bd0)
  xmax <- max(0.7, x9 + sign(x9) * port * bd9)

  inter <- 1
  n.add <- 32
  yleft <- rep(NA, n.add)
  yright <- rep(NA, n.add)

  if (xmax > x9) {
    idmost <- which.max(d$x)
    b.over <- (d$v[idmost] - d$v[idmost - 1]) / (d$x[idmost] - d$x[idmost - 1])
    extra.over <- seq(d$x[idmost], xmax, by = (xmax - d$x[idmost]) / n.add)
    b.over <- b.over * 1.2
    for (i in 1:length(extra.over)) {
      yright[i] <- d$v[idmost] - b.over * (d$x[idmost] - extra.over[i])
    }
  }

  if (x0 > xmin) {
    idlst <- which.min(d$x)
    b.bef <- (d$v[idlst + 1] - d$v[idlst]) / (d$x[idlst + 1] - d$x[idlst])
    extra.bef <- seq(xmin, d$x[idlst], by = (d$x[idlst] - xmin) / n.add)
    b.bef <- b.bef * 1.2
    for (i in 1:(length(extra.bef))) {
      yleft[i] <- d$v[idlst] - b.bef * (d$x[idlst] - extra.bef[i])
    }
  }

  msp1 <- data.frame(x = c(extra.bef, d$x, extra.over), y = c(yleft, d$v, yright))
  dv1 <- msp1$y

  w1 <- comp.prof.ci.new(d = msp1, v = dv1, conf = conf)

  w$w1 <- w1
  inter <- 1
  d <- msp1
  d$v <- msp1$y

  if (figure == TRUE & pick >= .90) {
    halfchi <- 0.5 * qchisq(conf[1], 1)
    w$ymin <- w$w1$vmax - (halfchi * 2.7) + 0.2
    w$ymax <- w$w1$vmax + 0.8 * halfchi - 0.2

    xlim0 <- -(w$w1$ci2 + (w$w1$ci_length) / 4)
    xlim1 <- -(w$w1$ci1 - (w$w1$ci_length) / 4)

    if (inter == 1) {
      co <- "green"
    } else if (inter == 2) {
      co <- "blue"
    }

    plot(-d$x, d$v,
         xlim = c(xlim0[1] - 0.02, xlim1[1] + 0.02),
         ylim = c(w$w1$vmax - (halfchi * 2.7), w$w1$vmax + 0.8 * halfchi),
         type = "l", xlab = "xi_Hosking in GEV",
         col = co, lwd = 2,
         ylab = "Profile Log-likelihood", main = c("Profile CI for xi"))
    ma <- w$w1$vmax
    abline(h = ma, col = 4)
    abline(h = ma - 0.5 * qchisq(conf[1], 1), col = 4)
    abline(v = c(-w$w1$xmax), col = 2, lty = 2, lwd = 2)
    text(x = c(-mle[3]), y = w$ymax, labels = "mle")
    abline(v = c(-w$w1$ci1[1], -w$w1$ci2[1]), col = 2, lty = 2)
    text(x = -w$w1$ci1[1], y = w$ymax,
         labels = paste("", round(-w$w1$ci1[1], 2), sep = ""), col = 2)
    text(x = -w$w1$ci2[1], y = w$ymax,
         labels = paste("", round(-w$w1$ci2[1], 2), sep = ""), col = 2)

    lm <- lmoms(data)
    lme <- pargev(lm)$para[3]
    abline(v = lme, lty = 5, lwd = 2, col = "purple")
    text(x = lme, y = w$ymin, labels = "lme")
  }

  return(w)
}


#' @keywords internal
comp.prof.ci.new <- function(d, v, conf = NULL) {
  w <- list()
  numcf <- length(conf)
  d$v <- v
  w$vmax <- max(d$v)
  nmsp <- length(d$v)

  w$nllh <- w$vmax
  idmax <- which.max(d$v)
  w$xmax <- d$x[idmax]

  ci1 <- rep(NA, numcf)
  ci2 <- rep(NA, numcf)
  ci_length <- rep(NA, numcf)

  for (i in 1:numcf) {
    chisq <- w$vmax - 0.5 * qchisq(conf[i], 1)
    diff <- d$v - chisq

    ci1[i] <- d$x[which.min(abs(diff[1:idmax]))]
    ci2[i] <- d$x[(idmax + 1):nmsp][which.min(abs(diff[(idmax + 1):nmsp]))]

    ci_min <- min(ci1[i], ci2[i])
    ci_max <- max(ci1[i], ci2[i])
    ci1[i] <- ci_min
    ci2[i] <- ci_max

    ci_length[i] <- ci2[i] - ci1[i]
  }

  w$ci1 <- ci1
  w$ci2 <- ci2
  w$ci_length <- ci_length

  return(w)
}
