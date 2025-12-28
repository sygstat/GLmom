# GLmomentEst 1.0.0

* Initial CRAN release.
* Stationary GEV parameter estimation using generalized L-moments (`glme.gev()`).
* Non-stationary GEV11 model estimation (`glme.gev11()`) where location
  (mu) and scale (sigma) parameters vary linearly with time.
  - `para.glme`: Proposed GLME estimates
  - `para.jkss`: Estimates based on Shin et al. (2025) JKSS paper
* Multiple penalty function options: `"beta"` (default), `"norm"`, `"ms"`
  (Martins-Stedinger), `"park"`, `"cannon"`, `"cd"` (Coles-Dixon), and
  `"no"` (no penalty).
* Flexible hyperparameter specification via `pen.choice` or direct parameters
  (`p`, `c1`, `c2` for beta; `mu`, `std` for normal penalty).
* Included datasets: `streamflow` and `PhliuAgromet`.
* References:
  - Shin et al. (2025) arXiv:2512.20385 (GLME methodology)
  - Shin et al. (2025) J. Korean Stat. Soc. 54, 947-970 (non-stationary L-moment)
