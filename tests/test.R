# GLmom package - basic tests
library(GLmom)

set.seed(777)

# --- 1. Data loading ---
data(streamflow)
data(Trehafod)
data(PhliuAgromet)
data(bangkok)
data(haenam)
data(glanteifi)

stopifnot(is.data.frame(streamflow), nrow(streamflow) == 50)
stopifnot(is.data.frame(Trehafod), nrow(Trehafod) == 53)
stopifnot(is.data.frame(PhliuAgromet), nrow(PhliuAgromet) == 40)
stopifnot(is.data.frame(bangkok), nrow(bangkok) == 58)
stopifnot(is.data.frame(haenam), nrow(haenam) == 52)
stopifnot(is.data.frame(glanteifi), nrow(glanteifi) == 65)
stopifnot(all(c("id", "river", "location", "year", "flow") %in% names(glanteifi)))

# --- 2. Stationary GEV (glme.gev) ---
res1 <- glme.gev(streamflow$r1)

stopifnot(is.list(res1))
stopifnot("para.glme" %in% names(res1))
stopifnot("para.lme" %in% names(res1))
stopifnot("nllh.glme" %in% names(res1))
stopifnot("convergence" %in% names(res1))
stopifnot(length(res1$para.glme) == 3)
stopifnot(length(res1$para.lme) == 3)
stopifnot(is.numeric(res1$para.glme))
stopifnot(res1$para.glme[2] > 0)  # sigma > 0
# v1.x compatibility fields still present
stopifnot(all(c("pen", "c1_c2", "covinv.lmom", "lcovdet") %in% names(res1)))

# --- 3. Penalty functions ---
for (pen in c("beta", "norm", "ms", "park", "cannon", "cd", "no")) {
  r <- glme.gev(streamflow$r1, pen = pen)
  stopifnot(is.numeric(r$para.glme), length(r$para.glme) == 3)
}

# All preset choices run
for (ch in 1:6) {
  r <- glme.gev(streamflow$r1, ntry = 3, pen = "beta", pen.choice = ch)
  stopifnot(is.numeric(r$para.glme))
}
for (ch in 1:4) {
  r <- glme.gev(streamflow$r1, ntry = 3, pen = "norm", pen.choice = ch)
  stopifnot(is.numeric(r$para.glme))
}

# v1.x-style call with explicit hyperparameters (pen.choice = NULL)
r_old <- glme.gev(streamflow$r1, ntry = 3, pen = "beta", pen.choice = NULL,
                  p = 6, c1 = 10, c2 = 5)
stopifnot(is.numeric(r_old$para.glme))

# pk.beta / pk.beta.stnary are the same implementation
pb1 <- pk.beta(para = c(100, 20, -0.2), lme.center = c(100, 20, -0.15),
               p = 6, c1 = 3, c2 = 1)
pb2 <- pk.beta.stnary(para = c(100, 20, -0.2), lme.center = c(100, 20, -0.15),
                      p = 6, c1 = 3, c2 = 1)
stopifnot(identical(pb1, pb2), pb1$pk.one > 0)

# --- 4. Non-stationary GEV11 (glme.gev11) ---
res2 <- glme.gev11(glanteifi$flow, ntry = 5)

stopifnot(is.list(res2))
stopifnot("para.glme" %in% names(res2))
stopifnot("convergence" %in% names(res2))
stopifnot(length(res2$para.glme) == 5)
# v1.x compatibility fields
stopifnot(all(c("strup.org", "para.wls", "para.gado", "lme.sta") %in% names(res2)))
stopifnot(length(res2$para.wls) == 5, length(res2$lme.sta) == 3)

# pen="no" returns para.lme and equals lme.gev11() (same RNG state)
set.seed(123)
res2n <- glme.gev11(glanteifi$flow, ntry = 5, pen = "no")
stopifnot("para.lme" %in% names(res2n))
set.seed(123)
lg <- lme.gev11(glanteifi$flow, ntry = 5)
# same L-moment solution up to solver/RNG tolerance (lmrob consumes RNG)
stopifnot(max(abs(res2n$para.glme - lg$lme.gev11)) < 1e-3)

# Regression baseline: LME of Glanteifi matches Shin et al. (2025), Table 4
# (LME row: 175.48, 0.6, 3.78, 0.005, -0.316)
stopifnot(abs(lg$lme.gev11[1] - 175.48) < 0.05)
stopifnot(abs(lg$lme.gev11[2] - 0.60)   < 0.01)
stopifnot(abs(lg$lme.gev11[3] - 3.78)   < 0.01)
stopifnot(abs(lg$lme.gev11[4] - 0.005)  < 0.001)
stopifnot(abs(lg$lme.gev11[5] - (-0.316)) < 0.005)

# opt.choose = "gof" also runs
res2g <- glme.gev11(glanteifi$flow, ntry = 5, opt.choose = "gof")
stopifnot(length(res2g$para.glme) == 5)

# --- 5. Component methods (new in v2.0.0) ---
s <- strup.gev11(glanteifi$flow)
stopifnot(all(c("strup.para", "strup.mdfy") %in% names(s)))
stopifnot(length(s$strup.mdfy) == 5)

gn <- GN16.gev11(glanteifi$flow)
stopifnot(all(c("para.gado.org", "para.gado.mdfy") %in% names(gn)))
stopifnot(length(gn$para.gado.org) == 5)

rr <- ran.gev_all(30, para = c(100, 0.5, 3, 0.005, -0.2), model = "gev11")
stopifnot(is.numeric(rr), length(rr) == 30, all(is.finite(rr)))

# --- 6. Compatibility wrappers (nsgev, gado.prop_11) ---
res3 <- nsgev(Trehafod$r1, ntry = 5)

stopifnot(is.list(res3))
stopifnot("para.prop" %in% names(res3))
stopifnot(length(res3$para.prop) == 5)

res3b <- suppressWarnings(gado.prop_11(Trehafod$r1, ntry = 5))
stopifnot(all(c("para.prop", "para.gado", "para.wls", "strup.org", "lme.sta")
              %in% names(res3b)))

# --- 7. Model averaging (ma.gev) ---
res4 <- ma.gev(streamflow$r1, quant = c(0.99, 0.995), numk = 7)

stopifnot(is.list(res4))
stopifnot("qua.ma" %in% names(res4))
stopifnot("qua.mle" %in% names(res4))
stopifnot("qua.lme" %in% names(res4))
stopifnot("w.ma" %in% names(res4))
stopifnot("pick_xi" %in% names(res4))
stopifnot(length(res4$qua.ma) == 2)
stopifnot(all(res4$qua.ma > 0))

# --- 8. Auxiliary functions ---
# pargev.kfix
lmom <- lmomco::lmoms(streamflow$r1)
pk <- pargev.kfix(lmom, kfix = -0.1)
stopifnot(is.list(pk), "para" %in% names(pk))

# quagev.NS (returns one quantile per time point)
q <- quagev.NS(f = 0.99, para = res2$para.glme, nsample = 65, model = "gev11")
stopifnot(is.numeric(q), length(q) == 65, all(q > 0))

# init.glme: new interface and v1.x-style named xdat=
i1 <- init.glme(streamflow$r1, ntry = 4)
i2 <- init.glme(xdat = streamflow$r1, ntry = 4)
stopifnot(is.matrix(i1), nrow(i1) == 4, is.matrix(i2), nrow(i2) == 4)
i3 <- init.glme(glanteifi$flow, ntry = 4, model = "gev11",
                pretheta = s$strup.mdfy)
stopifnot(is.matrix(i3), nrow(i3) == 4)

message("All tests passed.")
