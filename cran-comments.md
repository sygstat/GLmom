## GLmom 2.0.1

This is a correction release to 2.0.0 (published 2026-07-11). It fixes
two implementation errors identified by the methodology authors:

* The constant term of the GLME objective now uses log(det(V))/2 (the
  correct multivariate-normal constant) instead of log(det(V)). Parameter
  estimates are unchanged (the term is parameter-free); only reported
  criterion values shift.
* The fixed literature penalties (`pen = "ms"`, `"park"`, `"cannon"`) are
  now evaluated on the shape-parameter interval (-0.5, 0.5), matching
  their definitions in the literature (Martins & Stedinger, 2000), instead
  of (-1, 0.5). Estimates obtained with these three penalty options
  change; this affects reported results, which is why we submit the
  correction promptly rather than waiting the customary interval between
  updates.

Both changes are documented in NEWS.md.

## Test environments

* local: macOS (Apple silicon, arm64), R 4.1.3 — 0 errors, 0 warnings,
  2 notes (see below)
* win-builder, R-devel (2026-07-23 r90295 ucrt): Status OK —
  0 errors, 0 warnings, 0 notes
* macOS builder not used: the package is pure R and was fully checked
  locally on macOS (arm64)

## R CMD check results

0 errors | 0 warnings | 2 notes (local only; win-builder was clean)

1. **Days since last update**
   2.0.1 follows 2.0.0 after a short interval because it corrects the
   estimation results of three penalty options (see above); we considered
   it better to fix published results promptly.

2. **(possibly) invalid URL / DOI (status 403)**
   The flagged DOIs (10.1029/1999WR900330, AGU;
   10.1111/j.2517-6161.1990.tb01775.x, Wiley) are valid; the 403
   responses come from publishers rejecting automated requests.
   The links resolve correctly in a browser.

(A local "unable to verify current time" note also appears on the build
machine; it is an environment artifact unrelated to the package.)

## Downstream dependencies

There are no reverse dependencies on CRAN
(checked with `tools::package_dependencies("GLmom", reverse = TRUE)`).
