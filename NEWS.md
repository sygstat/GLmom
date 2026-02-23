# GLmom 1.3.1

## Bug Fixes

* Replaced all `\dontrun{}` with unwrapped examples or `\donttest{}` per CRAN policy.
* Wrapped `Trehafod` and `glme.gev11()` examples in `\donttest{}` (> 5 sec on Debian).
* Fixed invalid URL for UK National River Flow Archive Peak Flow Dataset.
* Single-quoted acronyms and proper names in DESCRIPTION to avoid spell-check notes.

# GLmom 1.3.0

## Breaking Changes

* **BREAKING**: `glme.gev()` output names changed for consistency:
  - `lme` → `para.lme`
  - `glme` → `para.glme`
  - `nllh.pref` → `nllh.glme`
  - `covinv` → `covinv.lmom`
* **BREAKING**: `ma.gev()` output names changed for consistency:
  - `zp.ma` → `qua.ma`
  - `zp.bma` → `qua.bma`
  - `fin.se.ma` → `fixw.se.ma`
  - `adj.se.ma` → `ranw.se.ma`
  - `numk_ma` and `numk_bma` → `run.numk`
  - `pick_xi_ma` and `pick_xi_bma` → `pick_xi`
  - `remle1` → `para.remle1` (in return and internal use)
  - `remle2` → `para.remle2` (in return and internal use)

## Internal Function Renames

* Renamed internal functions for consistency (these are not exported but may affect code using `:::` accessor):
  - `init.glme()` → `init.gevmax()`
  - `new_pf_norm()` → `pk.norm.stnary()` (with backward-compatible alias)
  - `gev.rl.delta_new()` → `gev.rl.delta()`
  - `lme.boots.new()` → `lme.boots()`
  - `cand.xi.new.paper()` → `cand.xi()`
  - `weight.com.new()` → `weight.com()`
  - `cov.interp.new()` → `cov.interp()`
  - `gev.profxi.mdfy.paper()` → `gev.profxi.mdfy()`
  - `comp.prof.ci.new()` → `comp.prof.ci()`
  - `gev1.CD()` → `mle.gev.CD()`
  - `gev.remle()` → `remle.gev()`
  - `ginit.max()` → `init.gevmax()`
  - `pargev.xifix.ma()` → `pargev.xifix()`
* Penalty functions in `set.prior()` now use `pk.beta.stnary()` from glme.gev.R

# GLmom 1.2.0

## Breaking Changes

* **BREAKING**: `glme.gev11()` output `para.jkss` renamed to `para.lme` for consistency.
  Users accessing `result$para.jkss` should update their code to use `result$para.lme`.
* **BREAKING**: `glme.gev11()` and `gado.prop_11()` output `strup.final` renamed to `para.wls`.
  Users accessing `result$strup.final` should update their code to use `result$para.wls`.
* **BREAKING**: `glme.gev11()` no longer returns `strup.sta` in its output.

## New Features

* Enhanced `glme.gev11()` with new parameters:
  - `glme.pre = "wls"`: Pre-estimation method selection ("wls" or "gado")
  - `choose = "gof"`: Model selection criterion ("gof" for goodness-of-fit, "nllh" for negative log-likelihood)
  - `pen.choice = 6`: Default penalty hyperparameter choice changed from NULL to 6
* New `quagev.NS()` function for calculating quantiles from non-stationary GEV models
  - Supports GEV11, GEV10, GEV20, and stationary GEV00 models
  - Returns time-varying quantiles as vector or matrix
* Enhanced `ma.gev()` with new estimation options:
  - `CD = TRUE`: Coles-Dixon penalized MLE for shape parameter regularization
  - `remle = TRUE`: Restricted MLE with mean/median constraints
  - Returns `mle.CD`, `qua.CD`, `remle1`, `remle2`, `qua.remle1`, `qua.remle2`
  - Returns `quant` in output for convenience
  - BMA outputs now include `bma.se.between` and `bma.se.within`
* New diagnostic plotting functions for MAGEV:
  - `magev.ksensplot()`: K sensitivity analysis to select optimal number of submodels
  - `magev.qqplot()`: 2x2 Q-Q diagnostic plot comparing MLE, LME, surrogate, and REMLE
  - `magev.rlplot()`: Return level plot with 95% confidence intervals
* Added `bangkok` dataset: Annual maximum daily rainfall from Bangkok, Thailand
* Added `haenam` dataset: Annual maximum daily rainfall from Haenam, South Korea

## Bug Fixes

* Improved handling of single quantile (`numq = 1`) in `ma.gev()`

# GLmom 1.1.0

* Added Model Averaging GEV estimation (`ma.gev()`) for high quantile estimation.
  - Combines MLE and L-moment estimates through weighted model averaging
  - Multiple weighting schemes: `like`, `gLd`, `med`, `cvt` and variants
  - Optional Bayesian model averaging (`bma=TRUE`) with normal/beta priors
  - Returns model-averaged quantiles (`qua.ma`) with standard errors
* New dependencies: `ismev`, `Rsolnp`, `zoo`.
* Reference: Shin et al. (2026) SERA, 40(2), 47 (MAGEV method)

# GLmom 1.0.0

* Initial CRAN release.
* Stationary GEV parameter estimation using generalized L-moments (`glme.gev()`).
* Non-stationary GEV11 model estimation (`glme.gev11()`) where location
  (mu) and scale (sigma) parameters vary linearly with time.
  - `para.glme`: Proposed GLME estimates
  - `para.lme`: L-moment based estimates for non-stationary model
* Compatibility functions for Shin et al. (2025, J. Korean Stat. Soc.):
  - `nsgev()`: Simple interface for L-moment based non-stationary estimation
  - `gado.prop_11()`: Comprehensive estimation with multiple methods
* Multiple penalty function options: `"beta"` (default), `"norm"`, `"ms"`
  (Martins-Stedinger), `"park"`, `"cannon"`, `"cd"` (Coles-Dixon), and
  `"no"` (no penalty).
* Flexible hyperparameter specification via `pen.choice` or direct parameters
  (`p`, `c1`, `c2` for beta; `mu`, `std` for normal penalty).
* Included datasets: `streamflow`, `PhliuAgromet`, and `Trehafod`.
* References:
  - Shin et al. (2025a) arXiv:2512.20385 (GLME method)
  - Shin et al. (2025b) J. Korean Stat. Soc. 54:947-970 (Non-stationary L-moment)
