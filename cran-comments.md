## R CMD check results

0 errors | 0 warnings | 3 notes

## Test environments

* Docker Debian (aarch64-unknown-linux-gnu), R 4.3.2

## NOTEs

1. **New submission**
   ```
   Maintainer: 'Yonggwan Shin <syg.stat@etri.re.kr>'
   New submission
   ```
   This is the first submission of this package to CRAN.

2. **Time verification**
   ```
   unable to verify current time
   ```
   This is a Docker environment issue, not a package problem.

3. **Pandoc not installed**
   ```
   Files 'README.md' or 'NEWS.md' cannot be checked without 'pandoc' being installed.
   ```
   This is a build environment issue. README.md renders correctly.

## Downstream dependencies

There are currently no downstream dependencies as this is a new package.

## Package description

This package provides generalized L-moments estimation methods for the
generalized extreme value (GEV) distribution:

- **Stationary GEV**: `glme.gev()` estimates GEV parameters with various
  penalty functions for shape parameter regularization.
- **Non-stationary GEV11**: `glme.gev11()` estimates time-varying GEV
  parameters where location mu(t) and scale sigma(t) change linearly with time.
- **Model Averaging (MAGEV)**: `ma.gev()` combines MLE and L-moment estimates
  through weighted model averaging for robust high quantile estimation.
- **Penalty functions**: beta (adaptive), normal, Martins-Stedinger, Park,
  Cannon, Coles-Dixon, and no penalty options.
- **Weighting schemes** (for MAGEV): likelihood-based (AIC), generalized
  L-moment distance, median-based, and conventional AIC.

The methodology is described in:
- Shin et al. (2025a) arXiv:2512.20385 (GLME method)
- Shin et al. (2025b) J. Korean Stat. Soc. 54:947-970 (Non-stationary L-moment)
- Shin et al. (2026) SERA, 40(2), 47 (MAGEV method)
