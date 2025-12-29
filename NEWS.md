# GLmomentEst 1.0.0

* Initial CRAN release.
* Stationary GEV parameter estimation using generalized L-moments (`glme.gev()`).
* Non-stationary GEV11 model estimation (`glme.gev11()`) where location
  (mu) and scale (sigma) parameters vary linearly with time.
  - `para.glme`: Proposed GLME estimates
  - `para.jkss`: L-moment based estimates for non-stationary model
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
