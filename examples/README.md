# GLmom Package Examples

Reproducible examples for the GLmom R package demonstrating Generalized L-Moment Estimation for Extreme Value Distributions.

## Quick Start

```r
# Install package (if not already installed)
remotes::install_github("sygstat/GLmom")

# Run all examples
setwd("path/to/GLmom/examples")
source("00_run_all_examples.R")

# Or run individual examples
source("01_stationary_gev.R")
```

## Example Files

| File | Description | Main Function |
|------|-------------|---------------|
| `01_stationary_gev.R` | Stationary GEV estimation with penalty functions | `glme.gev()` |
| `02_nonstationary_gev11.R` | Non-stationary GEV with time-varying parameters | `glme.gev11()` |
| `03_model_averaging.R` | Model averaging for high quantile estimation | `ma.gev()` |
| `04_compatibility_functions.R` | Shin et al. (2025b) methodology wrappers | `nsgev()`, `gado.prop_11()` |
| `05_magev_diagnostics.R` | MAGEV diagnostic plots and new datasets | `magev.ksensplot()`, `magev.qqplot()`, `magev.rlplot()` |
| `00_run_all_examples.R` | Run all examples sequentially | - |

## Example 1: Stationary GEV

Demonstrates parameter estimation for stationary GEV distribution:
- Default adaptive beta penalty
- Comparison of 7 penalty functions (beta, norm, ms, park, cannon, cd, no)
- Custom hyperparameter specification
- Quantile/return level estimation

```r
library(GLmom)
data(streamflow)
result <- glme.gev(streamflow$r1, pen = "beta")
result$glme  # (mu, sigma, xi)
```

## Example 2: Non-stationary GEV11

Demonstrates time-varying GEV estimation:
- Model: mu(t) = mu0 + mu1*t, sigma(t) = exp(sigma0 + sigma1*t), xi constant
- Multi-start optimization with `ntry`
- Trend interpretation
- Time-varying return levels

```r
library(GLmom)
data(Trehafod)
result <- glme.gev11(Trehafod$r1, ntry = 10)
result$para.glme  # (mu0, mu1, sigma0, sigma1, xi)
result$para.lme   # Pure L-moment estimates (no penalty)
```

**Note (v1.2.0)**: The output `para.jkss` has been renamed to `para.lme` for consistency. Update existing code accordingly.

## Example 3: Model Averaging (MAGEV)

Demonstrates model averaging for robust high quantile estimation:
- Multiple weighting schemes (like, gLd, med, cvt)
- Bootstrap standard errors
- Bayesian Model Averaging option
- Comparison of MLE, LME, and MA estimates

**New in v1.2.0:**
- `CD = TRUE`: Coles-Dixon penalized MLE for shape parameter regularization
- `remle = TRUE`: Restricted MLE with mean/median constraints
- Returns `quant` in output for convenience
- BMA outputs include `bma.se.between` and `bma.se.within`

```r
library(GLmom)
data(streamflow)

# Basic model averaging
result <- ma.gev(streamflow$r1, quant = c(0.99), weight = "like1", B = 200)
result$zp.ma  # Model-averaged quantile

# With CD and REMLE options (new in v1.2.0)
result <- ma.gev(streamflow$r1, quant = c(0.98, 0.99, 0.995),
                 weight = "like1", B = 100, CD = TRUE, remle = TRUE)
result$qua.CD      # CD-penalized MLE quantiles
result$qua.remle1  # REMLE (mean constraint) quantiles
result$qua.remle2  # REMLE (median constraint) quantiles
```

## Example 4: Compatibility Functions

Demonstrates simplified interfaces for Shin et al. (2025b) methodology:
- `nsgev()`: Returns only proposed L-moment estimates
- `gado.prop_11()`: Returns estimates from WLS, GN16, and proposed methods

```r
library(GLmom)
data(Trehafod)
result <- nsgev(Trehafod$r1, ntry = 10)
result$para.prop  # Proposed estimates
```

## Example 5: MAGEV Diagnostic Visualization (New in v1.2.0)

Demonstrates diagnostic plotting functions for MAGEV analysis:

### K Sensitivity Analysis
`magev.ksensplot()` helps determine the optimal number of submodels K:
```r
data(bangkok)
optimal_k <- magev.ksensplot(data = bangkok[,1], mink = 4, maxk = 20,
                              quant = c(0.99, 0.995))
```

### Q-Q Diagnostic Plot
`magev.qqplot()` creates a 2x2 panel comparing MLE, LME, Surrogate, and REMLE:
```r
qq <- c(seq(0.01, 0.99, by = 0.01), 0.995, 0.999)
zx <- ma.gev(bangkok[,1], quant = qq, weight = 'like1',
             numk = 9, remle = TRUE)
magev.qqplot(data = bangkok[,1], zx = zx)
```

### Return Level Plot
`magev.rlplot()` displays fitted return levels with 95% confidence intervals:
```r
ff <- c(seq(0.01, 0.09, by = 0.01), seq(0.1, 0.9, by = 0.1),
        0.93, 0.95, 0.98, 0.99, 0.995, 0.999)
zx <- ma.gev(bangkok[,1], quant = ff, weight = 'like1',
             numk = 9, varcom = TRUE)
magev.rlplot(par = zx$surr$par, se.vec = zx$adj.se.ma, data = bangkok[,1])
```

## Datasets

| Dataset | Description | n | Source |
|---------|-------------|---|--------|
| `streamflow` | Annual maximum streamflow | 50 | Hydrological data |
| `Trehafod` | River flow, Wales, UK (1968-2020) | 53 | UK National River Flow Archive |
| `PhliuAgromet` | Thai meteorological data | - | Phliu Agrometeorological Station |
| `bangkok` | Annual max daily rainfall, Bangkok, Thailand | - | Thai Meteorological Department (TMD) |
| `haenam` | Annual max daily rainfall, Haenam, South Korea | - | Korea Meteorological Administration (KMA) |

**Note**: `bangkok` and `haenam` datasets are new in v1.2.0.

## New Features in v1.2.0

### Breaking Changes
- `glme.gev11()` output `para.jkss` renamed to `para.lme` for consistency

### New Features
- **ma.gev() enhancements:**
  - `CD = TRUE`: Coles-Dixon penalized MLE
  - `remle = TRUE`: Restricted MLE with mean/median constraints
  - Returns `quant` in output
  - BMA outputs include `bma.se.between` and `bma.se.within`

- **New diagnostic plotting functions:**
  - `magev.ksensplot()`: K sensitivity analysis
  - `magev.qqplot()`: Q-Q diagnostic plot (2x2 panel)
  - `magev.rlplot()`: Return level plot with 95% CI

- **New datasets:**
  - `bangkok`: Annual maximum daily rainfall from Bangkok, Thailand
  - `haenam`: Annual maximum daily rainfall from Haenam, South Korea

## References

- Shin, Y., Shin, Y., Park, J., & Park, J.-S. (2025a). Generalized method of L-moment estimation for stationary and nonstationary extreme value models. *arXiv:2512.20385*

- Shin, Y., Shin, Y., & Park, J.-S. (2025b). Building nonstationary extreme value model using L-moments. *Journal of the Korean Statistical Society*, 54, 947-970.

- Shin, Y., Shin, Y., & Park, J.-S. (2025c). Model averaging with mixed criteria for estimating high quantiles of extreme values. *arXiv:2505.21417*

## Requirements

- R >= 3.5
- GLmom package
- Dependencies: lmomco, nleqslv, robustbase, ismev, Rsolnp, zoo

## License

GPL (>= 3)

## Contact

- Email: syg.stat@etri.re.kr
- GitHub: https://github.com/sygstat/GLmom
