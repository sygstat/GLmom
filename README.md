# GLmom

An R package for **Generalized L-moments Estimation** of the Generalized Extreme Value (GEV) Distribution.

## Overview

This package provides three main estimation approaches for extreme value analysis:

1. **GLME (Generalized L-Moment Estimation)**: Combines L-moments (Hosking, 1990) with penalty functions to regularize the shape parameter, providing more stable estimates especially for small samples (Shin et al., 2025a).
2. **NS L-moment based estimation**: Pure L-moment equations for non-stationary models without penalty (Shin et al., 2025b).
3. **MAGEV (Model Averaging GEV)**: Combines MLE and L-moment estimates through weighted model averaging for robust high quantile estimation (Shin et al., 2026).

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
- **Coles-Dixon penalized MLE** (`mle.gev.CD`): MLE with exponential penalty on shape parameter (Coles & Dixon, 1999), available via `ma.gev(..., CD=TRUE)`
- **Restricted MLE** (`remle.gev`): Constrained MLE matching sample mean (or median) and L-scale, available via `ma.gev(..., remle=TRUE)`
- **Shin et al. (2025b) compatibility** (`nsgev`, `gado.prop_11`): Pure L-moment based estimation
- **Multiple penalty functions**: `beta`, `norm`, `ms`, `park`, `cannon`, `cd`, `no`
- **Multiple weighting schemes** (MAGEV): `like`, `gLd`, `med`, `cvt`
- **MAGEV diagnostic plots**: `magev.ksensplot`, `magev.qqplot`, `magev.rlplot`
- **Example datasets**: `streamflow`, `PhliuAgromet`, `Trehafod`, `bangkok`, `haenam`

## Installation

```r
# Install from GitHub
remotes::install_github("sygstat/GLmom")
```

## Quick Start

```r
library(GLmom)

# Stationary GEV
data(streamflow)
result <- glme.gev(streamflow$r1)
result$para.glme  # GLME estimates: (mu, sigma, xi)

# Non-stationary GEV11
data(Trehafod)
result <- glme.gev11(Trehafod$r1)
result$para.glme  # (mu0, mu1, sigma0, sigma1, xi)
```

## Examples

### 1. Stationary GEV Estimation

```r
library(GLmom)
data(streamflow)
x <- streamflow$r1

# Default: beta penalty (adaptive)
result <- glme.gev(x)
result$para.glme  # GLME estimates
#> [1]  55.5346   8.7206  -0.4508
result$para.lme   # Traditional L-moment estimates
#>       xi    alpha    kappa
#>  55.4021   9.5229  -0.4075

# Compare different penalty functions
glme.gev(x, pen = "beta")$para.glme[3]    # xi = -0.4508
glme.gev(x, pen = "ms")$para.glme[3]      # xi = -0.2855 (Martins-Stedinger)
glme.gev(x, pen = "park")$para.glme[3]    # xi = -0.3474
glme.gev(x, pen = "no")$para.glme[3]      # xi = -0.4075 (no penalty = L-moment)
```

### 2. Non-stationary GEV11 Estimation

```r
library(GLmom)
data(Trehafod)
x <- Trehafod$r1  # 53 years of river flow data

# Estimate with GLME (default: beta penalty)
result <- glme.gev11(x, ntry = 10)

# GLME estimates (with penalty)
result$para.glme
#>        mu0        mu1     sigma0     sigma1         xi
#>  79.312877   1.243031   2.788088   0.013611  -0.087837

# L-moment estimates (no penalty) - Shin et al. (2025b)
result$para.lme
#>        mu0        mu1     sigma0     sigma1         xi
#>  79.330671   1.243031   2.792161   0.013611  -0.076654

# Interpretation:
# - mu0 = 79.31: baseline location at t=0
# - mu1 = 1.24: location increases ~1.2 units per year
# - sigma0, sigma1: log-scale parameters
# - xi = -0.09: bounded upper tail (Weibull-type)
```

### 3. Shin et al. (2025b) Compatibility Functions

For users following the methodology in Shin et al. (2025b):

```r
library(GLmom)
data(Trehafod)

# Simple interface - returns proposed L-moment estimates
result1 <- nsgev(Trehafod$r1, ntry = 10)
result1$para.prop
#>        mu0        mu1     sigma0     sigma1         xi
#>  79.330671   1.243031   2.792161   0.013611  -0.076654

# Comprehensive output - multiple estimation methods
result2 <- gado.prop_11(Trehafod$r1, ntry = 10)

# Compare methods:
result2$para.prop    # Proposed L-moment method
result2$para.gado    # GN16 method
result2$para.wls     # Weighted Least Squares
result2$lme.sta      # Stationary L-moments (mu, sigma, xi)
```

### 4. Comparing Penalty Functions

```r
library(GLmom)
data(Trehafod)

# All penalty options for non-stationary model
penalties <- c("beta", "norm", "ms", "park", "cannon", "cd", "no")

results <- sapply(penalties, function(p) {
  r <- glme.gev11(Trehafod$r1, ntry = 10, pen = p)
  r$para.glme[5]  # shape parameter xi
})

print(round(results, 4))
#>  beta.xi norm.xi   ms.xi park.xi cannon.xi  cd.xi   no.xi
#>  -0.0878 -0.0904 -0.0895 -0.0692   -0.0901 -0.0667 -0.0767
```

### 5. Custom Hyperparameters

```r
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

```r
library(GLmom)
data(streamflow)
x <- streamflow$r1

# Model averaging with likelihood weights (default)
result <- ma.gev(x, quant = c(0.95, 0.99, 0.995), weight = 'like1', B = 200)

# Compare estimates
result$qua.mle    # MLE quantiles
#> [1]  107.28  147.02  167.35
result$qua.lme    # L-moment quantiles
#> [1]  110.43  184.35  234.27
result$qua.ma     # Model-averaged quantiles (recommended)
#> [1]  109.55  187.10  241.77

# Standard errors
result$fixw.se.ma # SE under fixed weights
result$ranw.se.ma # SE under random weights

# Using generalized L-moment distance weights
result2 <- ma.gev(x, quant = c(0.99), weight = 'gLd')
print(result2$w.ma)  # Model weights across K submodels

# Using Bayesian Model Averaging
result3 <- ma.gev(x, quant = c(0.99), bma = TRUE, pen = "norm")
print(result3$qua.bma)
```

### 7. Coles-Dixon Penalized MLE and Restricted MLE

The `ma.gev()` function provides optional access to two additional estimation
methods via `CD` and `remle` arguments.

**Coles-Dixon penalized MLE** (`mle.gev.CD`) applies an exponential penalty on the shape parameter to prevent extreme negative estimates, following Coles & Dixon (1999).

**Restricted MLE** (`remle.gev`) maximizes the GEV likelihood subject to moment constraints. Two stages are computed:

- **Stage 1** (`remle1`): Single constraint — the GEV mean (L1) matches the sample mean
- **Stage 2** (`remle2`): Dual constraints — both L1 and L2 (L-scale) match the sample L-moments

```r
library(GLmom)
data(streamflow)
x <- streamflow$r1

# CD-MLE: Coles-Dixon penalized MLE
result_cd <- ma.gev(x, quant = c(0.99, 0.995), weight = "like1",
                    B = 50, CD = TRUE)
result_cd$mle.CD   # Coles-Dixon penalized MLE parameters (mu, sigma, xi)
result_cd$qua.CD   # Quantiles from CD-MLE

# REMLE: Restricted MLE with mean constraint
result_re <- ma.gev(x, quant = c(0.99, 0.995), weight = "like1",
                    B = 50, remle = TRUE)
result_re$remle1      # REMLE stage 1 parameters (mean constraint)
result_re$qua.remle1  # Quantiles from REMLE stage 1
result_re$remle2      # REMLE stage 2 parameters (mean + L-scale constraints)
result_re$qua.remle2  # Quantiles from REMLE stage 2

# Combined: all methods at once
result_all <- ma.gev(x, quant = c(0.98, 0.99, 0.995), weight = "like1",
                     B = 50, CD = TRUE, remle = TRUE)
result_all$qua.mle     # Standard MLE quantiles
result_all$qua.CD      # CD-penalized MLE quantiles
result_all$qua.remle1  # REMLE (mean) quantiles
result_all$qua.remle2  # REMLE (mean + L-scale) quantiles
result_all$qua.ma      # Model-averaged quantiles
```

## Datasets

| Dataset          | Description                        | n  | Variables |
| ---------------- | ---------------------------------- | -- | --------- |
| `streamflow`   | Annual maximum streamflow          | 50 | Year, r1  |
| `PhliuAgromet` | Meteorological data from Thailand  | -  | prec, ... |
| `Trehafod`     | River flow from Wales, UK          | 53 | Year, r1  |
| `bangkok`      | Annual max daily rainfall, Bangkok | 52 | rainfall  |
| `haenam`       | Annual max daily rainfall, Haenam  | 49 | rainfall  |

```r
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

| Function              | Description                        | Output                                                                          |
| --------------------- | ---------------------------------- | ------------------------------------------------------------------------------- |
| `glme.gev()`        | Stationary GEV estimation          | `para.glme`, `para.lme`, `nllh.glme`                                      |
| `glme.gev11()`      | Non-stationary GEV11               | `para.glme`, `para.lme`, `para.gado`, `para.wls`, ...                   |
| `ma.gev()`          | Model averaging for high quantiles | `qua.ma`, `qua.mle`, `qua.lme`, `w.ma`, `qua.CD`, `qua.remle1`, ... |
| `nsgev()`           | Simple L-moment interface          | `para.prop`, `precis`                                                       |
| `gado.prop_11()`    | Comprehensive L-moment             | `para.prop`, `para.gado`, `para.wls`, `lme.sta`                         |
| `quagev.NS()`       | NS GEV quantile function           | quantiles (vector/matrix)                                                       |
| `magev.ksensplot()` | K sensitivity plot                 | optimal K value                                                                 |
| `magev.qqplot()`    | Q-Q diagnostic plot                | (graphical)                                                                     |
| `magev.rlplot()`    | Return level plot                  | (graphical)                                                                     |

### Auxiliary Estimation Functions (via `ma.gev()`)

| Function         | Description                            | Activated by                | Output                                                     |
| ---------------- | -------------------------------------- | --------------------------- | ---------------------------------------------------------- |
| `mle.gev.CD()` | Coles-Dixon penalized MLE for GEV      | `ma.gev(..., CD=TRUE)`    | `$mle.CD`, `$qua.CD`                                   |
| `remle.gev()`  | Restricted MLE with moment constraints | `ma.gev(..., remle=TRUE)` | `$remle1`, `$qua.remle1`, `$remle2`, `$qua.remle2` |

- **`mle.gev.CD()`**: Maximizes the GEV log-likelihood with an exponential penalty on the shape parameter ξ (Coles & Dixon, 1999). The penalty is neutral for ξ ≥ 0 and increasingly penalizes extreme negative ξ values, preventing unrealistic bounded upper tails.
- **`remle.gev()`**: Maximizes the GEV log-likelihood subject to equality constraints. Stage 1 constrains the distribution mean (L1) to match the sample mean. Stage 2 additionally constrains L-scale (L2) to match the sample L-scale. Uses `Rsolnp::solnp()` for constrained optimization with multiple random starting points.

### Penalty Functions

| Penalty    | Description             | Parameters            | Reference                  |
| ---------- | ----------------------- | --------------------- | -------------------------- |
| `beta`   | Adaptive beta (default) | `p`, `c1`, `c2` | Shin et al. (2025a)        |
| `norm`   | Normal distribution     | `mu`, `std`       | -                          |
| `ms`     | Martins-Stedinger       | Beta(6,9) fixed       | Martins & Stedinger (2000) |
| `park`   | Park                    | Beta(2.5,2.5) fixed   | -                          |
| `cannon` | Cannon                  | Beta(2,3.3) fixed     | -                          |
| `cd`     | Coles-Dixon             | exponential           | Coles & Dixon (1999)       |
| `no`     | No penalty              | -                     | Pure L-moments             |

### Weighting Schemes (MAGEV)

| Weight                                | Description                   | Trim       |
| ------------------------------------- | ----------------------------- | ---------- |
| `like`, `like0`, `like1`        | Likelihood-based (AIC)        | 0, 0, 1    |
| `gLd`, `gLd0`, `gLd1`, `gLd2` | Generalized L-moment distance | 0, 0, 1, 2 |
| `med`, `med1`, `med2`           | Median + L-moment distance    | 0, 1, 2    |
| `cvt`                               | Conventional AIC              | -          |

### Helper Functions

| Function             | Description                              |
| -------------------- | ---------------------------------------- |
| `glme.like()`      | GLME likelihood function                 |
| `init.gevmax()`    | Parameter initialization                 |
| `pargev.kfix()`    | GEV with fixed shape                     |
| `MS_pk()`          | Martins-Stedinger penalty                |
| `pk.beta.stnary()` | Beta penalty function                    |
| `pk.norm.stnary()` | Normal penalty function                  |
| `cd.hos()`         | Coles-Dixon penalty for `mle.gev.CD()` |

## Authors

- **Yonggwan Shin**, Senior Researcher, Electronics and Telecommunications Research Institute, Korea ([syg.stat@etri.re.kr](mailto:syg.stat@etri.re.kr))
- **Yire Shin**, Ph.D, Chonnam National University, Korea([shinyire87@gmail.com](mailto:shinyire87@gmail.com))
- **Jihong Park**, Chonnam National University, Korea([jihong8090@gmail.com](mailto:jihong8090@gmail.com))
- **Jeong-Soo Park**, Professor, Chonnam National University, Korea([jspark@chonnam.ac.kr](mailto:jspark@chonnam.ac.kr))

## Citation

If you use this package, please cite:

> Shin, Y., Shin, Y., Park, J., & Park, J.-S. (2025). Generalized method of L-moment estimation for stationary and nonstationary extreme value models. *arXiv preprint* arXiv:2512.20385. https://doi.org/10.48550/arXiv.2512.20385
>
> Shin, Y., Shin, Y., & Park, J. S. (2026). Model averaging with mixed criteria for estimating high quantiles of extreme values: Application to heavy rainfall. *Stochastic Environmental Research and Risk Assessment*, 40(2), 47. https://doi.org/10.1007/s00477-025-03167-x

## References

- Shin, Y., Shin, Y., Park, J., & Park, J.-S. (2025a). Generalized method of L-moment estimation for stationary and nonstationary extreme value models. *arXiv preprint* arXiv:2512.20385.
- Shin, Y., Shin, Y., & Park, J.-S. (2025b). Building nonstationary extreme value model using L-moments. *Journal of the Korean Statistical Society*, 54, 947-970.
- Shin, Y., Shin, Y., & Park, J. S. (2026). Model averaging with mixed criteria for estimating high quantiles of extreme values: Application to heavy rainfall. *Stochastic Environmental Research and Risk Assessment*, 40(2), 47. https://doi.org/10.1007/s00477-025-03167-x
- Strupczewski, W. G., & Kaczmarek, Z. (2001). Non-stationary approach to at-site flood frequency modelling II. Weighted least squares estimation.  *Journal of Hydrology* ,  *248* (1-4), 143-151.
- Hosking, J.R.M. (1990). L-moments: Analysis and estimation of distributions using linear combinations of order statistics. *Journal of the Royal Statistical Society B*, 52, 105-124.
- Martins, E.S., & Stedinger, J.R. (2000). Generalized maximum-likelihood generalized extreme-value quantile estimators for hydrologic data. *Water Resources Research*, 36, 737-744.
- Coles, S., & Dixon, M. (1999). Likelihood-based inference for extreme value models. *Extremes*, 2, 5-23.

## License

GPL (>= 3)
