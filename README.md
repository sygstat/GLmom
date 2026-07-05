# GLmom

An R package for **Generalized L-moments Estimation** of the Generalized Extreme Value (GEV) Distribution.

## Overview

This package provides three main estimation approaches for extreme value analysis:

1. **GLME (Generalized L-Moment Estimation)**: Combines L-moments (Hosking, 1990) with penalty functions to regularize the shape parameter, providing more stable estimates especially for small samples (Shin et al., 2025a).
2. **NS L-moment based estimation**: Pure L-moment equations for non-stationary models without penalty (Shin et al., 2025b).
3. **MAGEV (Model Averaging GEV)**: Combines MLE and L-moment estimates through weighted model averaging for robust high quantile estimation (Shin et al., 2026).

### The GEV Distribution

The Generalized Extreme Value (GEV) distribution has CDF (Hosking's parameterization):

```
F(x) = exp{-[1 - xi*(x-mu)/sigma]^(1/xi)}
```

where `mu` is location, `sigma > 0` is scale, and `xi` is the shape parameter.
Note the sign convention: `xi < 0` corresponds to a heavy (Fréchet-type) upper
tail and `xi > 0` to a bounded (Weibull-type) upper tail. This is Hosking's
L-moment convention, used by `lmom`/`lmomco`; the shape parameter reported by
`evd`, `ismev`, and `extRemes` has the opposite sign.

### Non-stationary GEV11 Model

For time-varying extremes, the GEV11 model allows:

- **Location**: `mu(t) = mu0 + mu1 * t`
- **Scale**: `sigma(t) = exp(sigma0 + sigma1 * t)`
- **Shape**: `xi` (constant over time)

## Features

- **Stationary GEV estimation** (`glme.gev`): GLME with various penalty functions
- **Non-stationary GEV11 estimation** (`glme.gev11`): Time-varying parameters with GLME
- **Pure L-moment NS estimation** (`lme.gev11`, new in v2.0.0): Shin et al. (2025b) method (renames `gado.prop_11`)
- **Component NS methods** (new in v2.0.0): `strup.gev11` (Strupczewski WLS), `GN16.gev11` (quantile-based GN16)
- **NS random number generation** (`ran.gev_all`, new in v2.0.0) and quantiles (`quagev.NS`)
- **Model averaging estimation** (`ma.gev`): High quantile estimation with multiple weighting schemes
- **Coles-Dixon penalized MLE** (`mle.gev.CD`): MLE with exponential penalty on shape parameter (Coles & Dixon, 1999), available via `ma.gev(..., CD=TRUE)`
- **Restricted MLE** (`remle.gev`): Constrained MLE matching sample mean (or median) and L-scale, available via `ma.gev(..., remle=TRUE)`
- **Shin et al. (2025b) compatibility** (`nsgev`, `gado.prop_11`): Pure L-moment based estimation
- **Multiple penalty functions**: `beta`, `norm`, `ms`, `park`, `cannon`, `cd`, `no`
- **Multiple weighting schemes** (MAGEV): `like`, `gLd`, `med`, `cvt`
- **Classed fits with methods** (new in v2.0.0): `print()`, `summary()`, and `plot()` (Q-Q diagnostics) for the objects returned by `glme.gev`, `glme.gev11`, `lme.gev11`, and `ma.gev`
- **MAGEV diagnostic plots**: `magev.ksensplot`, `magev.qqplot`, `magev.rlplot`
- **Example datasets**: `PhliuAgromet`, `bangkok`, `haenam`

> **Note (v2.0.0):** the GLME penalty presets and several defaults were
> re-tuned following the revised paper (Shin et al., 2025a), so numerical
> results differ from v1.3.1. All v1.3.1 functions remain available.
> The `streamflow`, `Trehafod`, and `glanteifi` datasets were removed
> because their source (UK National River Flow Archive) does not permit
> redistribution; see `NEWS.md` for details.

## Installation

This README describes **version 2.0.0**. Until CRAN is updated (the CRAN
release is currently 1.3.1, which still ships the old datasets, presets,
and defaults), install 2.0.0 from GitHub:

```r
# Version 2.0.0 (this README)
remotes::install_github("sygstat/GLmom")
```

```r
# CRAN release (currently 1.3.1)
install.packages("GLmom")
```


## Quick Start

```r
library(GLmom)

# Stationary GEV
data(haenam)
result <- glme.gev(haenam$X1)
result            # classed object: print/summary/plot methods (v2.0.0)
result$para.glme  # fields remain directly accessible: (mu, sigma, xi)

# Non-stationary GEV11
data(PhliuAgromet)
result <- glme.gev11(PhliuAgromet$prec)
result$para.glme  # (mu0, mu1, sigma0, sigma1, xi)

# Pure L-moment method (Shin et al. 2025b)
result_lme <- lme.gev11(PhliuAgromet$prec)
result_lme$lme.gev11
```

## Examples

### 1. Stationary GEV Estimation

```r
library(GLmom)
data(haenam)
x <- haenam$X1  # 52 annual maxima of daily rainfall, Haenam, Korea

# Default: beta penalty (adaptive), preset pen.choice = 1
result <- glme.gev(x)
result$para.glme  # GLME estimates
#>       mu      sig       xi
#> 113.3014  37.3662  -0.3225
result$para.lme   # Traditional L-moment estimates
#>       mu      sig       xi
#> 113.4524  37.3533  -0.3104

# Compare different penalty functions
glme.gev(x, pen = "beta")$para.glme[3]    # xi ~= -0.323
glme.gev(x, pen = "ms")$para.glme[3]      # xi ~= -0.328 (Martins-Stedinger)
glme.gev(x, pen = "park")$para.glme[3]    # xi ~= -0.308
glme.gev(x, pen = "no")$para.glme[3]      # xi ~= -0.310 (no penalty = L-moment)
```

### 2. Non-stationary GEV11 Estimation

```r
library(GLmom)
data(PhliuAgromet)
x <- PhliuAgromet$prec  # 40 years of annual max daily precipitation, Thailand

# Estimate with GLME (default: beta penalty, pen.choice = 1)
set.seed(123)
result <- glme.gev11(x, ntry = 10)

# GLME estimates (with penalty)
result$para.glme
#>      mu0      mu1   sigma0   sigma1       xi
#> 126.9019   0.5700   2.7549   0.0408  -0.0356

# L-moment estimates (no penalty) - Shin et al. (2025b)
set.seed(123)
result_lme <- lme.gev11(x, ntry = 10)
result_lme$lme.gev11
#>      mu0      mu1   sigma0   sigma1       xi
#> 126.7720   0.5700   2.7591   0.0408  -0.0303

# Interpretation:
# - mu0 = 126.9: baseline location at t=0
# - mu1 = 0.57: location increases ~0.6 mm per year
# - sigma0, sigma1: log-scale parameters (scale also increases over time)
# - xi = -0.04: near-Gumbel upper tail (Hosking convention: xi < 0 is Frechet-type)

# Companion methods (new in v2.0.0)
strup.gev11(x)$strup.mdfy      # WLS (Strupczewski & Kaczmarek 2001)
GN16.gev11(x)$para.gado.org    # GN16 method
```

### 3. Shin et al. (2025b) Compatibility Functions

For users following the methodology in Shin et al. (2025b):

```r
library(GLmom)
data(PhliuAgromet)

# Recommended interface (v2.0.0): lme.gev11()
set.seed(123)
result0 <- lme.gev11(PhliuAgromet$prec, ntry = 10)
result0$lme.gev11

# Simple wrapper - returns proposed L-moment estimates
set.seed(123)
result1 <- nsgev(PhliuAgromet$prec, ntry = 10)
result1$para.prop
#>      mu0      mu1   sigma0   sigma1       xi
#> 126.7720   0.5700   2.7591   0.0408  -0.0303

# Comprehensive output - multiple estimation methods (deprecated wrapper)
result2 <- gado.prop_11(PhliuAgromet$prec, ntry = 10)

# Compare methods:
result2$para.prop    # Proposed L-moment method (= lme.gev11)
result2$para.gado    # GN16 method (= GN16.gev11()$para.gado.org)
result2$para.wls     # Weighted Least Squares (= strup.gev11()$strup.mdfy)
result2$lme.sta      # Stationary L-moments (mu, sigma, xi)
```

### 4. Comparing Penalty Functions

```r
library(GLmom)
data(PhliuAgromet)

# All penalty options for non-stationary model
penalties <- c("beta", "norm", "ms", "park", "cannon", "cd", "no")

results <- sapply(penalties, function(p) {
  set.seed(123)
  r <- glme.gev11(PhliuAgromet$prec, ntry = 10, pen = p)
  r$para.glme[5]  # shape parameter xi
})

print(round(results, 4))
#>  beta.xi norm.xi   ms.xi park.xi cannon.xi   cd.xi   no.xi
#>  -0.0356 -0.0451 -0.1289 -0.0416   -0.0577 -0.0232 -0.0303
```

### 5. Custom Hyperparameters

```r
library(GLmom)
data(haenam)

# Preset choices (1-6 for beta, 1-4 for norm; Table 1 / Eq. 15 of the paper)
glme.gev(haenam$X1, pen = "beta", pen.choice = 2)
glme.gev(haenam$X1, pen = "norm", pen.choice = 1)

# Beta penalty with custom hyperparameters (set pen.choice = NULL)
glme.gev(haenam$X1, pen = "beta", pen.choice = NULL, p = 6, c1 = 5, c2 = 2)

# Normal penalty with custom mean and std
glme.gev(haenam$X1, pen = "norm", pen.choice = NULL, mu = -0.5, std = 0.2)
```

### 6. Model Averaging for High Quantiles (MAGEV)

```r
library(GLmom)
data(haenam)
x <- haenam$X1

# Model averaging with likelihood weights (default)
set.seed(42)
result <- ma.gev(x, quant = c(0.95, 0.99, 0.995), weight = 'like1', B = 200)

# Compare estimates
result$qua.mle    # MLE quantiles
#> [1] 310.72 569.44 741.63
result$qua.lme    # L-moment quantiles
#> [1] 295.67 494.90 615.84
result$qua.ma     # Model-averaged quantiles (recommended)
#> [1] 296.38 517.61 659.83

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
data(haenam)
x <- haenam$X1

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

| Dataset          | Description                                        | n  | Variables |
| ---------------- | -------------------------------------------------- | -- | --------- |
| `PhliuAgromet` | Annual max daily precipitation, Phliu, Thailand (1984-2023, has trend) | 40 | year, prec, ... |
| `bangkok`      | Annual max daily rainfall, Bangkok                 | 58 | X1-X5     |
| `haenam`       | Annual max daily rainfall, Haenam, Korea           | 52 | year, X1  |

```r
# Load and explore datasets
data(PhliuAgromet)
head(PhliuAgromet$prec)
#> [1] 108.1 103.1 140.5 137.6 188.9 123.8
```

## Function Reference

### Main Functions

| Function              | Description                        | Output                                                                          |
| --------------------- | ---------------------------------- | ------------------------------------------------------------------------------- |
| `glme.gev()`        | Stationary GEV estimation          | `para.glme`, `para.lme`, `nllh.glme`, `convergence`                     |
| `glme.gev11()`      | Non-stationary GEV11 (GLME)        | `para.glme`, `nllh.glme`, `convergence`, `para.wls`, `para.gado`, ...  |
| `lme.gev11()`       | Pure L-moment NS estimation (v2.0.0) | `lme.gev11`, `precis`                                                     |
| `strup.gev11()`     | WLS estimation (Strupczewski)      | `strup.para`, `strup.mdfy`                                                  |
| `GN16.gev11()`      | GN16 estimation                    | `para.gado.org`, `para.gado.mdfy`                                           |
| `ran.gev_all()`     | NS GEV random number generation    | numeric vector                                                                  |
| `ma.gev()`          | Model averaging for high quantiles | `qua.ma`, `qua.mle`, `qua.lme`, `w.ma`, `qua.CD`, `qua.remle1`, ... |
| `nsgev()`           | Simple L-moment interface          | `para.prop`, `precis`                                                       |
| `gado.prop_11()`    | Comprehensive L-moment (deprecated) | `para.prop`, `para.gado`, `para.wls`, `lme.sta`                        |
| `quagev.NS()`       | NS GEV quantile function           | quantiles (vector/matrix)                                                       |
| `magev.ksensplot()` | K sensitivity plot                 | optimal K value                                                                 |
| `magev.qqplot()`    | Q-Q diagnostic plot                | (graphical)                                                                     |
| `magev.rlplot()`    | Return level plot                  | (graphical)                                                                     |
| `print()` / `summary()` / `plot()` | Methods for fitted objects (`glme`, `glme11`, `lme11`, `magev` classes) | console summary / Q-Q diagnostics |

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
| `beta`   | Adaptive beta (default) | `p`, `c0`, `c1`, `c2` | Shin et al. (2025a)  |
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
| `init.glme()`      | Multi-start initialization (gev00/gev11) |
| `init.gevmax()`    | Parameter initialization (v1.x)          |
| `pargev.kfix()`    | GEV with fixed shape                     |
| `MS_pk()`          | Martins-Stedinger penalty                |
| `pk.beta()`        | Unified adaptive beta penalty (v2.0.0; alias `pk.beta.stnary()`) |
| `pk.norm.stnary()` | Normal penalty function                  |
| `cd.hos()`         | Coles-Dixon penalty for `mle.gev.CD()` |

## Authors

- **Yonggwan Shin**, Senior Researcher, Electronics and Telecommunications Research Institute, Korea ([syg.stat@etri.re.kr](mailto:syg.stat@etri.re.kr))
- **Seokkap Ko**, Principal Researcher, Electronics and Telecommunications Research Institute, Korea ([softgear@etri.re.kr](mailto:softgear@etri.re.kr))
- **Jihong Park**, Chonnam National University, Korea ([jihong8090@gmail.com](mailto:jihong8090@gmail.com))
- **Yire Shin**, Ph.D, Chonnam National University, Korea ([shinyire87@gmail.com](mailto:shinyire87@gmail.com))
- **Jeong-Soo Park**, Professor, Chonnam National University, Korea ([jspark@chonnam.ac.kr](mailto:jspark@chonnam.ac.kr))

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
