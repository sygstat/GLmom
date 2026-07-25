## GLmom 2.0.1

This is a quick correction release to 2.0.0. It fixes two implementation
errors identified by the methodology authors:

* The constant term of the GLME objective now uses log(det(V))/2 (the
  correct multivariate-normal constant) instead of log(det(V)). Parameter
  estimates are unchanged (the term is parameter-free); only reported
  criterion values shift.
* The fixed literature penalties (`pen = "ms"`, `"park"`, `"cannon"`) are
  now evaluated on the shape-parameter interval (-0.5, 0.5), matching
  their definitions in the literature (Martins & Stedinger, 2000), instead
  of (-1, 0.5). Estimates obtained with these three penalty options
  change; this affects reported results, which is why we submit the
  correction promptly rather than waiting.

Both changes are documented in NEWS.md.

## Previous release: GLmom 2.0.0

This is a major update (1.3.1 -> 2.0.0). It integrates the revised GLME
methodology of Shin et al. (2026 revised version of arXiv:2512.20385):
new exported functions for the non-stationary GEV11 model (`lme.gev11()`,
`strup.gev11()`, `GN16.gev11()`, `ran.gev_all()`), a unified penalty layer
(`pk.beta()`), re-tuned penalty presets, and a substantially faster
penalized non-stationary fit. All changes are documented in NEWS.md.

### Removed datasets

The datasets `streamflow`, `Trehafod`, and `glanteifi` (the latter only
ever present in the development version) have been **removed** in this
release. All three originate from the UK National River Flow Archive,
whose data terms and conditions
(<https://eidc.ceh.ac.uk/licences/NRFA-Data-Terms-and-Conditions>)
do not permit making the data available for download or redistributing
them to third parties, so we cannot ship them in the package. The
documentation points users to the NRFA Peak Flow Dataset for direct
access. The removal is prominently documented in NEWS.md; the remaining
example datasets (`PhliuAgromet`, `bangkok`, `haenam`) cover all examples
and tests.

## Test environments

* local (2.0.1): macOS (Apple silicon, arm64), R 4.1.3 — 0 errors,
  0 warnings, 2 notes (see below)
* win-builder, R-devel (2026-07-23 r90295 ucrt): Status OK —
  0 errors, 0 warnings, 0 notes
* macOS builder: service unavailable at the 2.0.0 submission; the package
  is pure R and was fully checked locally on macOS (arm64)

## R CMD check results

0 errors | 0 warnings | 2 notes

0. **Days since last update**
   2.0.1 follows 2.0.0 after a short interval because it corrects the
   estimation results of three penalty options (see above); we considered
   it better to fix published results promptly.

1. **(possibly) invalid URL / DOI (status 403)**
   The flagged DOIs (e.g., 10.1029/1999WR900330, AGU;
   10.1111/j.2517-6161.1990.tb01775.x, Wiley) are valid; the 403
   responses come from publishers rejecting automated requests.
   The links resolve correctly in a browser.

2. **unable to verify current time**
   Local environment artifact (no network time check available on the
   build machine); not related to the package.

3. **Possibly misspelled words in DESCRIPTION** (win-builder):
   'Coles', 'Hosking', 'Stedinger' are author surnames cited in the
   Description field. (They were single-quoted in an earlier submission
   and unquoted at CRAN's request.)

## Downstream dependencies

There are no reverse dependencies on CRAN
(checked with `tools::package_dependencies("GLmom", reverse = TRUE)`).
The API changes in this major release therefore affect no other CRAN
packages. Backward compatibility for end users is retained where the
data license permits: every function exported in 1.3.1 remains exported
and callable (with `gado.prop_11()` formally deprecated in favor of
`lme.gev11()`).
