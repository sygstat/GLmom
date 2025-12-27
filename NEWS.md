# GLmomentEst 1.0.0

* Initial CRAN release.
* Stationary GEV parameter estimation using generalized L-moments (`glme.gev()`).
* Non-stationary GEV11 model estimation (`glme.gev11()`) where location
  (mu) and scale (sigma) parameters vary linearly with time.
* Multiple penalty function options: `"beta"` (default), `"norm"`, `"ms"`
  (Martins-Stedinger), `"park"`, `"cannon"`, `"cd"` (Coles-Dixon), and
  `"no"` (no penalty).
* Flexible hyperparameter specification via `pen.choice` or direct parameters
  (`p`, `c1`, `c2` for beta; `mu`, `std` for normal penalty).
* Included datasets: `streamflow` and `PhliuAgromet`.
