
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GLmomentEst

An R package for Generalized L-moments Estimation of the Generalized Extreme Value (GEV) Distribution.

## Features

- **Stationary GEV estimation** (`glme.gev`): Estimate GEV parameters using generalized L-moments with various penalty functions
- **Non-stationary GEV11 estimation** (`glme.gev11`): Estimate time-varying GEV parameters where location and scale change linearly with time
- **Multiple penalty functions**: `beta`, `norm`, `ms` (Martins-Stedinger), `park`, `cannon`, `cd`, and `no` (no penalty)

## Installation

Install the latest development version (on GitHub) via `{remotes}`:

``` r
remotes::install_github("sygstat/GL-momentEst")
```

## Getting started

### Stationary GEV estimation

``` r
library(GLmomentEst)

# Load example data
data("streamflow")

# Basic usage with default beta penalty
result <- glme.gev(streamflow$r1)
result$glme  # GLME estimates (location, scale, shape)
result$lme   # L-moment estimates for comparison

# With different penalty functions
glme.gev(streamflow$r1, pen = "norm", mu = -0.5, std = 0.2)
glme.gev(streamflow$r1, pen = "ms")      # Martins-Stedinger
glme.gev(streamflow$r1, pen = "park")    # Park
glme.gev(streamflow$r1, pen = "cannon")  # Cannon
glme.gev(streamflow$r1, pen = "cd")      # CD
glme.gev(streamflow$r1, pen = "no")      # No penalty

# Using preset hyperparameter choices
glme.gev(streamflow$r1, pen = "beta", pen.choice = 2)
glme.gev(streamflow$r1, pen = "norm", pen.choice = 3)
```

### Non-stationary GEV11 estimation

The GEV11 model assumes time-varying parameters:
- Location: `mu(t) = mu0 + mu1 * t`
- Scale: `sigma(t) = exp(sigma0 + sigma1 * t)`
- Shape: `xi` (constant)

``` r
library(GLmomentEst)

# Load example data
data("PhliuAgromet")

# Estimate non-stationary GEV11 model
result <- glme.gev11(PhliuAgromet$prec)

# Results
result$para.glme   # Proposed GLME estimates (mu0, mu1, sigma0, sigma1, xi)
result$para.jkss   # L-moment based estimates for non-stationary model
result$lme.sta     # Stationary L-moment estimates
result$para.gado   # GN16 original estimates
result$strup.final # Weighted least squares estimates
```

## Function Reference

### Main Functions

| Function | Description |
|----------|-------------|
| `glme.gev()` | Stationary GEV parameter estimation using GLME |
| `glme.gev11()` | Non-stationary GEV11 model estimation |
| `glme.like()` | Likelihood function for GLME optimization |
| `init.glme()` | Initialize parameters for optimization |

### Penalty Functions

| Penalty | Description | Hyperparameters |
|---------|-------------|-----------------|
| `beta` | Adaptive beta distribution (default) | `p`, `c1`, `c2` |
| `norm` | Normal distribution | `mu`, `std` |
| `ms` | Martins-Stedinger (Beta(6,9)) | - |
| `park` | Park (Beta(2.5,2.5)) | - |
| `cannon` | Cannon (Beta(2,3.3)) | - |
| `cd` | CD penalty | - |
| `no` | No penalty | - |

## Authors

- **Yonggwan Shin**, Ph.D, Research & Development Center, XRAI Inc., Gwangju 61186, Korea, [syg.stat@gmail.com](mailto:syg.stat@gmail.com)
- **Yire Shin**, Ph.D, Department of Statistics, Chonnam National University, Republic of Korea, [shinyire87@gmail.com](mailto:shinyire87@gmail.com)
- **Jihong Park**, Department of Mathematics and Statistics, Chonnam National University, Republic of Korea
- **Jeong-Soo Park**, Professor, Department of Statistics, Chonnam National University, Republic of Korea, [jspark@jnu.ac.kr](mailto:jspark@jnu.ac.kr)

## Citation

If you use this package, please cite:

> Shin, Y., Shin, Y., Park, J., & Park, J.-S. (2025). Generalized method of L-moment estimation for stationary and nonstationary extreme value models. *arXiv preprint* arXiv:2512.20385. https://doi.org/10.48550/arXiv.2512.20385

## License

GPL (>= 3)

-----
