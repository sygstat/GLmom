# GLmom

An R package for **Generalized L-moments Estimation** of the Generalized Extreme Value (GEV) Distribution.

## Overview

This package provides three main estimation approaches for extreme value analysis:

1. **GLME (Generalized L-Moment Estimation)**: Combines L-moments (Hosking, 1990) with penalty functions to regularize the shape parameter, providing more stable estimates especially for small samples (Shin et al., 2025a).

2. **NS L-moment based estimation**: Pure L-moment equations for non-stationary models without penalty (Shin et al., 2025b).

3. **MAGEV (Model Averaging GEV)**: Combines MLE and L-moment estimates through weighted model averaging for robust high quantile estimation (Shin et al., 2025c).

### The GEV Distribution

The Generalized Extreme Value (GEV) distribution has CDF:

```
F(x) = exp{-[1 + xi*(x-mu)/sigma]^(-1/xi)}
```

where `mu` is location, `sigma > 0` is scale, and `xi` is shape parameter.

### Non-stationary GEV11 Model

For time-varying extremes, the GEV11 model allows:
- **Location**: `mu(t) = mu0 + mu1 * t`
- **Scale**: `sigma(t) = exp(sigma0 + sigma1 * t)`
- **Shape**: `xi` (constant over time)

## Features

- **Stationary GEV estimation** (`glme.gev`): GLME with various penalty functions
- **Non-stationary GEV11 estimation** (`glme.gev11`): Time-varying parameters with GLME
- **Model averaging estimation** (`ma.gev`): High quantile estimation with multiple weighting schemes
- **Shin et al. (2025b) compatibility** (`nsgev`, `gado.prop_11`): Pure L-moment based estimation
- **Multiple penalty functions**: `beta`, `norm`, `ms`, `park`, `cannon`, `cd`, `no`
- **Multiple weighting schemes** (MAGEV): `like`, `gLd`, `med`, `cvt`
- **MAGEV diagnostic plots**: `magev.ksensplot`, `magev.qqplot`, `magev.rlplot`
- **Example datasets**: `streamflow`, `PhliuAgromet`, `Trehafod`, `bangkok`, `haenam`

## Installation

``` r
# Install from GitHub
remotes::install_github("sygstat/GLmom")
```

## Quick Start

``` r
library(GLmom)

# Stationary GEV
data(streamflow)
result <- glme.gev(streamflow$r1)
result$glme  # GLME estimates: (mu, sigma, xi)

# Non-stationary GEV11
data(Trehafod)
result <- glme.gev11(Trehafod$r1)
result$para.glme  # (mu0, mu1, sigma0, sigma1, xi)
```

## Examples

### 1. Stationary GEV Estimation

``` r
library(GLmom)
data(streamflow)
x <- streamflow$r1

# Default: beta penalty (adaptive)
result <- glme.gev(x)
result$glme  # GLME estimates
#> [1]  55.5348   8.7204  -0.4508
result$lme   # Traditional L-moment estimates
#> [1]  55.5286   8.8033  -0.4999

# Compare different penalty functions
glme.gev(x, pen = "beta")$glme[3]    # xi = -0.4508
glme.gev(x, pen = "ms")$glme[3]      # xi = -0.3678 (Martins-Stedinger)
glme.gev(x, pen = "park")$glme[3]    # xi = -0.4648
glme.gev(x, pen = "no")$glme[3]      # xi = -0.4999 (no penalty = L-moment)
```

### 2. Non-stationary GEV11 Estimation

``` r
library(GLmom)
data(Trehafod)
x <- Trehafod$r1  # 53 years of river flow data

# Estimate with GLME (default: beta penalty)
result <- glme.gev11(x, ntry = 10)

# GLME estimates (with penalty)
result$para.glme
#>        mu0        mu1     sigma0     sigma1         xi
#>  84.549596   1.027800   2.910336   0.008994  -0.089364

# L-moment estimates (no penalty) - Shin et al. (2025b)
result$para.lme
#>        mu0        mu1     sigma0     sigma1         xi
#>  84.556050   1.027800   2.915524   0.008994  -0.078462

# Interpretation:
# - mu0 = 84.55: baseline location at t=0
# - mu1 = 1.03: location increases ~1 unit per year
# - sigma0, sigma1: log-scale parameters
# - xi = -0.08: bounded upper tail (Weibull-type)
```

### 3. Shin et al. (2025b) Compatibility Functions

For users following the methodology in Shin et al. (2025b):

``` r
library(GLmom)
data(Trehafod)

# Simple interface - returns proposed L-moment estimates
result1 <- nsgev(Trehafod$r1, ntry = 10)
result1$para.prop
#>        mu0        mu1     sigma0     sigma1         xi
#>  84.556050   1.027800   2.915524   0.008994  -0.078462

# Comprehensive output - multiple estimation methods
result2 <- gado.prop_11(Trehafod$r1, ntry = 10)

# Compare methods:
result2$para.prop    # Proposed L-moment method
result2$para.gado    # GN16 method
result2$para.wls     # Weighted Least Squares
result2$lme.sta      # Stationary L-moments (mu, sigma, xi)
```

### 4. Comparing Penalty Functions

``` r
library(GLmom)
data(Trehafod)

# All penalty options for non-stationary model
penalties <- c("beta", "norm", "ms", "park", "cannon", "cd", "no")

results <- sapply(penalties, function(p) {
  r <- glme.gev11(Trehafod$r1, ntry = 10, pen = p)
  r$para.glme[5]  # shape parameter xi
})

print(round(results, 4))
#>   beta   norm     ms   park cannon     cd     no
#> -0.089 -0.091 -0.090 -0.071 -0.091 -0.069 -0.078
```

### 5. Custom Hyperparameters

``` r
library(GLmom)
data(streamflow)

# Beta penalty with custom hyperparameters
glme.gev(streamflow$r1, pen = "beta", p = 6, c1 = 20, c2 = 7)

# Or use preset choices (1-6 for beta, 1-4 for norm)
glme.gev(streamflow$r1, pen = "beta", pen.choice = 2)
glme.gev(streamflow$r1, pen = "norm", pen.choice = 1)

# Normal penalty with custom mean and std
glme.gev(streamflow$r1, pen = "norm", mu = -0.5, std = 0.2)
```

### 6. Model Averaging for High Quantiles (MAGEV)

``` r
library(GLmom)
data(streamflow)
x <- streamflow$r1

# Model averaging with likelihood weights (default)
result <- ma.gev(x, quant = c(0.95, 0.99, 0.995), weight = 'like1', B = 200)

# Compare estimates
result$qua.mle    # MLE quantiles
#> [1]  72.52  85.57  93.45
result$qua.lme    # L-moment quantiles
#> [1]  72.95  87.75  97.20
result$zp.ma      # Model-averaged quantiles (recommended)
#> [1]  72.78  86.52  94.89

# Standard errors
result$fin.se.ma  # SE under fixed weights
result$adj.se.ma  # SE under random weights

# Using generalized L-moment distance weights
result2 <- ma.gev(x, quant = c(0.99), weight = 'gLd')
print(result2$w.ma)  # Model weights across K submodels

# Using Bayesian Model Averaging
result3 <- ma.gev(x, quant = c(0.99), bma = TRUE, pen = "norm")
print(result3$zp.bma)
```

## Datasets

| Dataset | Description | n | Variables |
|---------|-------------|---|-----------|
| `streamflow` | Annual maximum streamflow | 50 | Year, r1 |
| `PhliuAgromet` | Meteorological data from Thailand | - | prec, ... |
| `Trehafod` | River flow from Wales, UK | 53 | Year, r1 |
| `bangkok` | Annual max daily rainfall, Bangkok | 52 | rainfall |
| `haenam` | Annual max daily rainfall, Haenam | 49 | rainfall |

``` r
# Load and explore datasets
data(Trehafod)
head(Trehafod)
#>   Year    r1
#> 1 1968 72.61
#> 2 1969 52.82
#> 3 1970 94.80
```

## Function Reference

### Main Functions

| Function | Description | Output |
|----------|-------------|--------|
| `glme.gev()` | Stationary GEV estimation | `glme`, `lme`, `nllh` |
| `glme.gev11()` | Non-stationary GEV11 | `para.glme`, `para.lme`, `para.gado`, `para.wls`, ... |
| `ma.gev()` | Model averaging for high quantiles | `zp.ma`, `qua.mle`, `qua.lme`, `w.ma`, `qua.CD`, `qua.remle1`, ... |
| `nsgev()` | Simple L-moment interface | `para.prop`, `precis` |
| `gado.prop_11()` | Comprehensive L-moment | `para.prop`, `para.gado`, `para.wls`, `lme.sta` |
| `quagev.NS()` | NS GEV quantile function | quantiles (vector/matrix) |
| `magev.ksensplot()` | K sensitivity plot | optimal K value |
| `magev.qqplot()` | Q-Q diagnostic plot | (graphical) |
| `magev.rlplot()` | Return level plot | (graphical) |

### Penalty Functions

| Penalty | Description | Parameters | Reference |
|---------|-------------|------------|-----------|
| `beta` | Adaptive beta (default) | `p`, `c1`, `c2` | Shin et al. (2025a) |
| `norm` | Normal distribution | `mu`, `std` | - |
| `ms` | Martins-Stedinger | Beta(6,9) fixed | Martins & Stedinger (2000) |
| `park` | Park | Beta(2.5,2.5) fixed | - |
| `cannon` | Cannon | Beta(2,3.3) fixed | - |
| `cd` | Coles-Dixon | exponential | Coles & Dixon (1999) |
| `no` | No penalty | - | Pure L-moments |

### Weighting Schemes (MAGEV)

| Weight | Description | Trim |
|--------|-------------|------|
| `like`, `like0`, `like1` | Likelihood-based (AIC) | 0, 0, 1 |
| `gLd`, `gLd0`, `gLd1`, `gLd2` | Generalized L-moment distance | 0, 0, 1, 2 |
| `med`, `med1`, `med2` | Median + L-moment distance | 0, 1, 2 |
| `cvt` | Conventional AIC | - |

### Helper Functions

| Function | Description |
|----------|-------------|
| `glme.like()` | GLME likelihood function |
| `init.glme()` | Parameter initialization |
| `pargev.kfix()` | GEV with fixed shape |
| `MS_pk()` | Martins-Stedinger penalty |
| `pk.beta.stnary()` | Beta penalty function |

## Authors

- **Yonggwan Shin**, Senior Researcher, Electronics and Telecommunications Research Institute, Korea ([syg.stat@etri.re.kr](mailto:syg.stat@etri.re.kr))
- **Yire Shin**, Ph.D, Chonnam National University, Korea
- **Jihong Park**, Chonnam National University, Korea
- **Jeong-Soo Park**, Professor, Chonnam National University, Korea

## Citation

If you use this package, please cite:

> Shin, Y., Shin, Y., Park, J., & Park, J.-S. (2025). Generalized method of L-moment estimation for stationary and nonstationary extreme value models. *arXiv preprint* arXiv:2512.20385. https://doi.org/10.48550/arXiv.2512.20385

## References

- Shin, Y., Shin, Y., Park, J., & Park, J.-S. (2025a). Generalized method of L-moment estimation for stationary and nonstationary extreme value models. *arXiv preprint* arXiv:2512.20385.
- Shin, Y., Shin, Y., & Park, J.-S. (2025b). Building nonstationary extreme value model using L-moments. *Journal of the Korean Statistical Society*, 54, 947-970.
- Shin, Y., Shin, Y., & Park, J.-S. (2025c). Model averaging with mixed criteria for estimating high quantiles of extreme values: Application to heavy rainfall. *arXiv preprint* arXiv:2505.21417.
- Hosking, J.R.M. (1990). L-moments: Analysis and estimation of distributions using linear combinations of order statistics. *Journal of the Royal Statistical Society B*, 52, 105-124.
- Martins, E.S., & Stedinger, J.R. (2000). Generalized maximum-likelihood generalized extreme-value quantile estimators for hydrologic data. *Water Resources Research*, 36, 737-744.
- Coles, S., & Dixon, M. (1999). Likelihood-based inference for extreme value models. *Extremes*, 2, 5-23.

## License

GPL (>= 3)
