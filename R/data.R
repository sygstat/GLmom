#' Phliu Agrometeorological Station Data
#'
#' Annual maximum daily precipitation and station metadata from the Phliu
#' Agrometeorological Station, Chanthaburi province, Thailand (1984-2023).
#' A Mann-Kendall test indicates a significant increasing trend in the
#' annual maxima (tau = 0.235, p = 0.033), so this dataset is used for
#' demonstrating non-stationary (GEV11) estimation; it is the heavy
#' rainfall application of Shin et al. (2025).
#'
#' @format A data frame with 40 rows and 9 columns:
#' \describe{
#'   \item{Station.ID}{Station identifier (character)}
#'   \item{year}{Year of observation (numeric, 1984-2023)}
#'   \item{prec}{Annual maximum daily precipitation in mm (numeric)}
#'   \item{Name}{Station name (character)}
#'   \item{zone}{Climate zone code (character)}
#'   \item{latitude}{Station latitude in degrees (numeric)}
#'   \item{longitude}{Station longitude in degrees (numeric)}
#'   \item{Starting.year}{Record start year (integer)}
#'   \item{Ending.year}{Record end year (numeric)}
#' }
#'
#' @source Phliu Agrometeorological Station, Thailand.
#' 
#' @references 
#' Shin, Y., Shin, Y., Park, J., & Park, J. S. (2025). Generalized method of L-moment estimation
#' for stationary and nonstationary extreme value models. arXiv preprint arXiv:2512.20385.
#' \doi{10.48550/arXiv.2512.20385}
#'
#' @examples
#' data(PhliuAgromet)
#' head(PhliuAgromet)
#'
#' \donttest{
#' # Fit non-stationary GEV11 model by GLME
#' result <- glme.gev11(PhliuAgromet$prec, ntry = 5)
#' print(result$para.glme)
#' }
#'
"PhliuAgromet"

#' Bangkok Maximum Rainfall Data
#'
#' Annual maximum daily rainfall data from Bangkok, Thailand.
#' This dataset is used for demonstrating model averaging methods
#' for high quantile estimation in extreme value analysis.
#'
#' @format A data frame with 58 rows and 5 columns:
#' \describe{
#'   \item{X1}{Annual maximum daily rainfall in mm (numeric)}
#'   \item{X2}{2nd largest annual daily rainfall in mm (numeric)}
#'   \item{X3}{3rd largest annual daily rainfall in mm (numeric)}
#'   \item{X4}{4th largest annual daily rainfall in mm (numeric)}
#'   \item{X5}{5th largest annual daily rainfall in mm (numeric)}
#' }
#'
#' @source Thai Meteorological Department (TMD; \url{https://www.tmd.go.th})
#'
#' @references
#' Shin, Y., Shin, Y., & Park, J. S. (2026). Model averaging with mixed criteria
#' for estimating high quantiles of extreme values: Application to heavy rainfall.
#' \emph{Stochastic Environmental Research and Risk Assessment}, 40(2), 47.
#' \doi{10.1007/s00477-025-03167-x}
#'
#' @examples
#' data(bangkok)
#' head(bangkok)
#'
#' # Estimate high quantiles using model averaging
#' result <- ma.gev(bangkok$X1, quant = c(0.99, 0.995))
#' print(result$qua.ma)
#'
"bangkok"

#' Haenam Maximum Rainfall Data
#'
#' Annual maximum daily rainfall data from Haenam, South Korea.
#' This dataset is used for demonstrating model averaging methods
#' for high quantile estimation in extreme value analysis.
#'
#' @format A data frame with 52 rows and 2 columns:
#' \describe{
#'   \item{year}{Year of observation (integer, 1971-2022)}
#'   \item{X1}{Annual maximum daily rainfall in mm (numeric)}
#' }
#'
#' @source Korea Meteorological Administration (KMA; \url{https://www.weather.go.kr})
#'
#' @references
#' Shin, Y., Shin, Y., & Park, J. S. (2026). Model averaging with mixed criteria
#' for estimating high quantiles of extreme values: Application to heavy rainfall.
#' \emph{Stochastic Environmental Research and Risk Assessment}, 40(2), 47.
#' \doi{10.1007/s00477-025-03167-x}
#'
#' @examples
#' data(haenam)
#' head(haenam)
#'
#' # Estimate high quantiles using model averaging
#' result <- ma.gev(haenam$X1, quant = c(0.98, 0.99, 0.995))
#' print(result$qua.ma)
#'
"haenam"
