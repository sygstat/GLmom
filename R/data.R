#' Streamflow Data
#'
#' Annual maximum streamflow measurements for extreme value analysis.
#'
#' @format A data frame with 50 rows and 2 columns:
#' \describe{
#'   \item{Year}{Year of observation (character)}
#'   \item{r1}{Annual maximum streamflow value (numeric)}
#' }
#'
#' @source Hydrological data for extreme value analysis.
#'
#' @examples
#' data(streamflow)
#' head(streamflow)
#'
"streamflow"

#' Phliu Agrometeorological Station Data
#'
#' Climate or meteorological data from the Phliu Agrometeorological Station
#' for extreme value analysis.
#'
#' @format A data frame containing meteorological measurements.
#'
#' @source Phliu Agrometeorological Station.
#'
#' @examples
#' data(PhliuAgromet)
#' head(PhliuAgromet)
#'
"PhliuAgromet"

#' Trehafod River Flow Data
#'
#' Annual maximum river flow data from the Trehafod gauging station in Wales, UK.
#' This dataset is commonly used for demonstrating non-stationary extreme value
#' analysis methods.
#'
#' @format A data frame with 53 rows and 2 columns:
#' \describe{
#'   \item{Year}{Year of observation (1968-2020)}
#'   \item{r1}{Annual maximum river flow in cubic meters per second (m^3/s)}
#' }
#'
#' @source UK National River Flow Archive.
#'
#' @references
#' Shin, Y., Shin, Y. & Park, J.-S. (2025). Building nonstationary extreme value
#' model using L-moments. Journal of the Korean Statistical Society, 54, 947-970.
#' \doi{10.1007/s42952-025-00325-3}
#'
#' @examples
#' data(Trehafod)
#' head(Trehafod)
#'
#' # Fit non-stationary GEV11 model
#' \donttest{
#' result <- glme.gev11(Trehafod$r1, ntry = 5)
#' print(result$para.glme)
#' }
#'
"Trehafod"

#' Bangkok Maximum Rainfall Data
#'
#' Annual maximum daily rainfall data from Bangkok, Thailand.
#' This dataset is used for demonstrating model averaging methods
#' for high quantile estimation in extreme value analysis.
#'
#' @format A data frame containing annual maximum daily rainfall values.
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
#' \donttest{
#' result <- ma.gev(bangkok[,1], quant = c(0.99, 0.995))
#' print(result$qua.ma)
#' }
#'
"bangkok"

#' Haenam Maximum Rainfall Data
#'
#' Annual maximum daily rainfall data from Haenam, South Korea.
#' This dataset is used for demonstrating model averaging methods
#' for high quantile estimation in extreme value analysis.
#'
#' @format A data frame containing annual maximum daily rainfall values.
#'
#' @source Korea Meteorological Administration (KMA; \url{https://www.kma.go.kr})
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
#' \donttest{
#' result <- ma.gev(haenam[,1], quant = c(0.98, 0.99, 0.995))
#' print(result$qua.ma)
#' }
#'
"haenam"
