## R CMD check results

0 errors | 0 warnings | 3 notes

## Test environments

* local macOS Sequoia 15.6.1 (aarch64-apple-darwin20), R 4.1.3

## NOTEs

1. **New submission**
   ```
   Maintainer: 'Yonggwan Shin <syg.stat@gmail.com>'
   New submission
   ```
   This is the first submission of this package to CRAN.

2. **DOI check**
   ```
   Found the following (possibly) invalid DOIs:
     DOI: 10.1111/j.2517-6161.1990.tb01775.x
       From: DESCRIPTION
       Status: Forbidden
       Message: 403
   ```
   The DOI for Hosking (1990) returns a 403 Forbidden status during automated
   checking, but the DOI is valid and resolves correctly in web browsers.
   This is a known issue with some publisher DOI servers blocking automated requests.

3. **Time verification**
   ```
   unable to verify current time
   ```
   This is a local network/environment issue, not a package problem.

## Downstream dependencies

There are currently no downstream dependencies as this is a new package.

## Package description

This package provides generalized L-moments estimation methods for the
generalized extreme value (GEV) distribution:

- **Stationary GEV**: `glme.gev()` estimates GEV parameters with various
  penalty functions for shape parameter regularization.
- **Non-stationary GEV11**: `glme.gev11()` estimates time-varying GEV
  parameters where location μ(t) and scale σ(t) change linearly with time.
- **Penalty functions**: beta (adaptive), normal, Martins-Stedinger, Park,
  Cannon, Coles-Dixon, and no penalty options.

The methodology is described in:
- Shin et al. (2025) arXiv:2512.20385
- Shin et al. (2025) J. Korean Stat. Soc. 54:947-970
