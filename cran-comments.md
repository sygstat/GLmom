## Resubmission

This is a resubmission. In this version I have:

* Replaced all `\dontrun{}` with unwrapped examples or `\donttest{}`.
* Wrapped `Trehafod` and `glme.gev11()` examples in `\donttest{}`
  as they exceed 5 sec on Debian.
* Fixed the invalid URL for the UK National River Flow Archive
  (`https://nrfa.ceh.ac.uk/peak-flow-dataset` changed to
  `https://nrfa.ceh.ac.uk/data/peak-flow-dataset`).
* Single-quoted acronyms and proper names in DESCRIPTION
  (e.g., 'GEV', 'GLME', 'Hosking', 'Stedinger') to address
  the misspelled words NOTE.

## R CMD check results

0 errors | 0 warnings | 1 note

## Test environments

* macOS Tahoe 26.2 (aarch64-apple-darwin20), R 4.5.2

## NOTEs

1. **New submission**
   ```
   Maintainer: 'Yonggwan Shin <syg.stat@etri.re.kr>'
   New submission
   ```
   This is the first submission of this package to CRAN.

## Downstream dependencies

There are currently no downstream dependencies as this is a new package.
