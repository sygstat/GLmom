# GLmom 2.0.0

This release integrates the revised GLME methodology of Shin et al.
(2026 revised version of arXiv:2512.20385) for both
the stationary GEV and the non-stationary GEV11 models. All functions and
datasets of v1.3.1 remain available and callable, but several defaults,
hyperparameter presets, and output fields changed, so numerical results
differ from v1.3.1.

## New Features

* New exported functions for the non-stationary GEV11 model:
  - `lme.gev11()` - pure L-moment estimation (renames `gado.prop_11()`;
    Shin et al. 2025b, JKSS).
  - `strup.gev11()` - weighted least squares estimation
    (Strupczewski & Kaczmarek 2001) with the modified final specification.
  - `GN16.gev11()` - quantile-based GN16 estimation.
  - `ran.gev_all()` - random number generation from stationary and
    non-stationary GEV models (gev, gev01, gev10, gev11, gev20).
* New exported penalty function `pk.beta()`, the unified data-adaptive beta
  penalty now shared by the stationary and non-stationary methods.
  `pk.beta.stnary()` is kept as an alias of the new implementation.
* The `PhliuAgromet` documentation now notes the significant increasing
  trend in the annual maxima (Mann-Kendall tau = 0.235, p = 0.033), which
  makes it the recommended example data for the non-stationary methods.
* Fitted objects are now classed (`"glme"`, `"glme11"`, `"lme11"`,
  `"magev"`) with `print()`, `summary()`, and `plot()` methods; `plot()`
  draws quantile-quantile diagnostics (on the standard Gumbel scale for
  the GEV11 fits). The objects remain plain lists, so all documented
  fields stay accessible; a `data` field storing the input series was
  added to support the diagnostics.
* `glme.gev()` and `glme.gev11()` gained arguments `c0` (beta penalty
  support half-width), `q` (fixed beta shape), and `show` (verbose);
  `glme.gev()` also gained `method`, `maxit`, `abstol` for `optim()`.

## Breaking Changes

* **The `streamflow`, `Trehafod`, and `glanteifi` datasets were removed.**
  All three originate from the UK National River Flow Archive (NRFA), whose
  data terms and conditions do not permit making the data available for
  download or redistributing them to third parties, so they cannot be
  shipped in a CRAN package. Users can obtain the underlying series
  directly from the NRFA Peak Flow Dataset
  (\<https://nrfa.ceh.ac.uk/peak-flow-dataset\>). The remaining example
  datasets are `PhliuAgromet`, `bangkok`, and `haenam`; all examples,
  tests, and documentation now use these.
* **Penalty hyperparameter presets were re-tuned** to match Table 1 and
  Eq. (15) of the revised paper:
  - beta `pen.choice` 1-6: (p, c1, c2) = (6,3,1), (6,5,2), (6,7,3),
    (2,3,0.5), (2,5,1), (2,7,1.5)  [was (6,10,5), (6,20,7), (6,30,9),
    (2,10,5), (2,20,7), (2,30,9)].
  - norm `pen.choice` 1-4: (mu, std) = (-0.5,0.25), (-0.5,0.15),
    (-0.6,0.25), (-0.6,0.15)  [std was 0.2/0.1].
* **Default arguments changed**: `glme.gev(ntry=5, pen.choice=1, c1=3, c2=1)`
  (was `ntry=10, pen.choice=NULL, c1=10, c2=5`);
  `glme.gev11(ntry=5, opt.choose="nllh", pen.choice=1, std=0.2, c1=5, c2=2)`
  (was `ntry=10, opt.choose="gof", pen.choice=NULL, std=0.3, c1=10, c2=5`).
* **`glme.gev11()` was restructured** following the revised methodology:
  the penalized estimation now starts from the WLS pre-estimate of
  `strup.gev11()` and uses the fixed asymptotic covariance matrix of the
  sample L-moments (Eq. 26 of the paper, rescaled by 300/n) instead of the
  bootstrap covariance. This is substantially faster.
  - Output field `para.lme` is now returned only when `pen="no"`
    (use `lme.gev11()` to obtain the pure L-moment estimates otherwise).
  - Output fields `strup.org`, `para.wls`, `para.gado`, `lme.sta` are
    still returned for compatibility.
  - New output fields: `convergence`, `pen_pen.choice`, `c0_c1_c2`
    (in addition to `pen`, `p_q`, `c1_c2`, `mu_std`).
  - Argument `glme.pre` is deprecated and ignored (warning);
    `init.rob` is still honored (passed to `strup.gev11()`).
* **`pk.beta.stnary()` behavior changed** (now an alias of `pk.beta()`):
  the adaptive support is (max(-1, xi-c0), min(0.3, xi+c0)) with default
  `c0=0.35` (was 0.3), the q-adaptation triggers at xi <= -0.05 (was 0),
  values outside the support return 1e-100 (was 1), and for fixed `q` the
  default support is (-1, 0.5) (was (-0.5, 0.5)). Consequently the fixed
  literature penalties ("ms", "park", "cannon") inside `glme.gev()` /
  `glme.gev11()` are now evaluated on (-1, 0.5), following the revised
  reference code.
* `init.glme()` is now the multi-model initializer (first argument `data`,
  new arguments `model`, `pretheta`); calls using the old named argument
  `xdat=` still work. The old behavior is available as `init.gevmax()`.
* `glme.gev()` output: parameter vectors are now named (mu, sig, xi) and
  `convergence`, `pen_pen.choice`, `c0_c1_c2` were added. The v1.x fields
  `pen`, `p_q`, `c1_c2`, `covinv.lmom`, `lcovdet` are still returned.

## Deprecated

* `gado.prop_11()` is deprecated in favor of `lme.gev11()`; it still
  returns the v1.x output fields (para.prop, para.gado, para.wls,
  strup.org, lme.sta), now reassembled from `lme.gev11()`,
  `GN16.gev11()`, and `strup.gev11()`.
* `nsgev()` is kept as a convenience wrapper around `lme.gev11()`.

## Internal Changes

* New internal functions: `check.penalty()`, `penalty.fun()` (unified
  penalty dispatcher), `boot.cov()`, `fun.lme.gev11()`, `find_max_beta.pk()`,
  `trans.gum01()` (vectorized).
* Removed internal functions (affects `:::` users only): `pk.beta.ns()`,
  `obj.lme.gev11()`, `gev11.GLD()`, `init.glme.gev11()`,
  `strup.glme.gev11()`, and the old `optim.glme.gev11()` (replaced by
  `opt.glme.gev11()`).
* MAGEV is untouched in this release: `set.prior()` (BMA prior) now uses a
  frozen internal copy of the v1.x beta preference function, so `ma.gev()`
  results are identical to v1.3.1.
* Robustness fixes in the multi-start loops: a failed `optim()`/`nleqslv()`
  try no longer aborts the remaining tries.
* `time.m.gev11()`: fixed the `checkmom=` typo (now `checklmom=FALSE`).
* Console output now uses `message()` or is gated behind `show=TRUE`.

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
