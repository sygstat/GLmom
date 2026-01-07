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
- Model: μ(t) = μ₀ + μ₁t, σ(t) = exp(σ₀ + σ₁t), ξ constant
- Multi-start optimization with `ntry`
- Trend interpretation
- Time-varying return levels

```r
library(GLmom)
data(Trehafod)
result <- glme.gev11(Trehafod$r1, ntry = 10)
result$para.glme  # (mu0, mu1, sigma0, sigma1, xi)
```

## Example 3: Model Averaging (MAGEV)

Demonstrates model averaging for robust high quantile estimation:
- Multiple weighting schemes (like, gLd, med, cvt)
- Bootstrap standard errors
- Bayesian Model Averaging option
- Comparison of MLE, LME, and MA estimates

```r
library(GLmom)
data(streamflow)
result <- ma.gev(streamflow$r1, quant = c(0.99), weight = "like1", B = 200)
result$zp.ma  # Model-averaged quantile
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

## Datasets

| Dataset | Description | n | Variables |
|---------|-------------|---|-----------|
| `streamflow` | Annual maximum streamflow | 50 | Year, r1 |
| `Trehafod` | River flow, Wales, UK (1968-2020) | 53 | Year, r1 |
| `PhliuAgromet` | Thai meteorological data | - | prec, ... |

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
