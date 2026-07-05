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
data(haenam)
result <- glme.gev(haenam$X1, pen = "beta")
result$para.glme  # (mu, sigma, xi)
```

## Example 2: Non-stationary GEV11

Demonstrates time-varying GEV estimation:
- Model: mu(t) = mu0 + mu1*t, sigma(t) = exp(sigma0 + sigma1*t), xi constant
- Multi-start optimization with `ntry`
- Trend interpretation
- Time-varying return levels

```r
library(GLmom)
data(PhliuAgromet)
result <- glme.gev11(PhliuAgromet$prec, ntry = 10)
result$para.glme  # (mu0, mu1, sigma0, sigma1, xi)
result$para.lme   # Pure L-moment estimates (no penalty)
```

## Example 3: Model Averaging (MAGEV)

Demonstrates model averaging for robust high quantile estimation:
- Multiple weighting schemes (like, gLd, med, cvt)
- Bootstrap standard errors
- Bayesian Model Averaging option
- Comparison of MLE, LME, and MA estimates

Additional features:
- `CD = TRUE`: Coles-Dixon penalized MLE for shape parameter regularization
- `remle = TRUE`: Restricted MLE with mean constraint (stage 1) and mean+L-scale constraints (stage 2)
- Returns `quant` in output for convenience
- BMA outputs include `bma.se.between` and `bma.se.within`

```r
library(GLmom)
data(haenam)

# Basic model averaging
result <- ma.gev(haenam$X1, quant = c(0.99), weight = "like1", B = 200)
result$qua.ma  # Model-averaged quantile

# With CD and REMLE options
result <- ma.gev(haenam$X1, quant = c(0.98, 0.99, 0.995),
                 weight = "like1", B = 100, CD = TRUE, remle = TRUE)
result$qua.CD      # CD-penalized MLE quantiles
result$qua.remle1  # REMLE (mean constraint) quantiles
result$qua.remle2  # REMLE (mean + L-scale constraints) quantiles
```

## Example 4: Compatibility Functions

Demonstrates simplified interfaces for Shin et al. (2025b) methodology:
- `nsgev()`: Returns only proposed L-moment estimates
- `gado.prop_11()`: Returns estimates from WLS, GN16, and proposed methods

```r
library(GLmom)
data(PhliuAgromet)
result <- nsgev(PhliuAgromet$prec, ntry = 10)
result$para.prop  # Proposed estimates
```

## Example 5: MAGEV Diagnostic Visualization

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
ff <- c(seq(0.01, 0.09, by = 0.01), 0.1, 0.2, 0.3, 0.4, 0.5,
        0.6, 0.7, 0.8, 0.9, 0.93, 0.95, 0.98, 0.99,
        0.993, 0.995, 0.998, 0.999)
zx <- ma.gev(bangkok[,1], quant = ff, weight = 'like1',
             numk = 9, varcom = TRUE)
magev.rlplot(par = zx$surr$par, se.vec = zx$ranw.se.ma, data = bangkok[,1])
```

## Datasets

| Dataset | Description | n | Source |
|---------|-------------|---|--------|

| `PhliuAgromet` | Thai meteorological data | 40 | Phliu Agrometeorological Station |
| `bangkok` | Annual max daily rainfall, Bangkok, Thailand | 58 | Thai Meteorological Department (TMD) |
| `haenam` | Annual max daily rainfall, Haenam, South Korea | 52 | Korea Meteorological Administration (KMA) |

## References

- Shin, Y., Shin, Y., Park, J., & Park, J.-S. (2025a). Generalized method of L-moment estimation for stationary and nonstationary extreme value models. *arXiv:2512.20385*

- Shin, Y., Shin, Y., & Park, J.-S. (2025b). Building nonstationary extreme value model using L-moments. *Journal of the Korean Statistical Society*, 54, 947-970.

- Shin, Y., Shin, Y., & Park, J. S. (2026). Model averaging with mixed criteria for estimating high quantiles of extreme values: Application to heavy rainfall. *Stochastic Environmental Research and Risk Assessment*, 40(2), 47. https://doi.org/10.1007/s00477-025-03167-x

## Fitted-object methods (v2.0.0)

The estimation functions return classed objects (`glme`, `glme11`, `lme11`,
`magev`) with `print()`, `summary()`, and `plot()` (Q-Q diagnostics) methods;
examples 01 and 02 demonstrate them.

## Requirements

- R >= 3.5
- GLmom package
- Dependencies: lmomco, nleqslv, robustbase, ismev, Rsolnp, zoo

## License

GPL (>= 3)

## Contact

- Email: syg.stat@etri.re.kr
- GitHub: https://github.com/sygstat/GLmom
