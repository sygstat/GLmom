# GLmom package - basic tests
library(GLmom)

# --- 1. Data loading ---
data(streamflow)
data(Trehafod)
data(PhliuAgromet)
data(bangkok)
data(haenam)

stopifnot(is.data.frame(streamflow), nrow(streamflow) == 50)
stopifnot(is.data.frame(Trehafod), nrow(Trehafod) == 53)
stopifnot(is.data.frame(PhliuAgromet), nrow(PhliuAgromet) == 40)
stopifnot(is.data.frame(bangkok), nrow(bangkok) == 58)
stopifnot(is.data.frame(haenam), nrow(haenam) == 52)

# --- 2. Stationary GEV (glme.gev) ---
res1 <- glme.gev(streamflow$r1)

stopifnot(is.list(res1))
stopifnot("para.glme" %in% names(res1))
stopifnot("para.lme" %in% names(res1))
stopifnot("nllh.glme" %in% names(res1))
stopifnot(length(res1$para.glme) == 3)
stopifnot(length(res1$para.lme) == 3)
stopifnot(is.numeric(res1$para.glme))
stopifnot(res1$para.glme[2] > 0)  # sigma > 0

# --- 3. Penalty functions ---
for (pen in c("beta", "norm", "ms", "park", "cannon", "cd", "no")) {
  r <- glme.gev(streamflow$r1, pen = pen)
  stopifnot(is.numeric(r$para.glme), length(r$para.glme) == 3)
}

# --- 4. Non-stationary GEV11 (glme.gev11) ---
res2 <- glme.gev11(Trehafod$r1, ntry = 5)

stopifnot(is.list(res2))
stopifnot("para.glme" %in% names(res2))
stopifnot("para.lme" %in% names(res2))
stopifnot("para.gado" %in% names(res2))
stopifnot(length(res2$para.glme) == 5)

# --- 5. Compatibility wrappers (nsgev, gado.prop_11) ---
res3 <- nsgev(Trehafod$r1, ntry = 5)

stopifnot(is.list(res3))
stopifnot("para.prop" %in% names(res3))
stopifnot(length(res3$para.prop) == 5)

# --- 6. Model averaging (ma.gev) ---
res4 <- ma.gev(streamflow$r1, quant = c(0.99, 0.995), numk = 7)

stopifnot(is.list(res4))
stopifnot("qua.ma" %in% names(res4))
stopifnot("qua.mle" %in% names(res4))
stopifnot("qua.lme" %in% names(res4))
stopifnot("w.ma" %in% names(res4))
stopifnot("pick_xi" %in% names(res4))
stopifnot(length(res4$qua.ma) == 2)
stopifnot(all(res4$qua.ma > 0))

# --- 7. Auxiliary functions ---
# pargev.kfix
lmom <- lmomco::lmoms(streamflow$r1)
pk <- pargev.kfix(lmom, kfix = -0.1)
stopifnot(is.list(pk), "para" %in% names(pk))

# quagev.NS (returns one quantile per time point)
q <- quagev.NS(f = 0.99, para = res2$para.glme, nsample = 53, model = "gev11")
stopifnot(is.numeric(q), length(q) == 53, all(q > 0))

message("All tests passed.")
