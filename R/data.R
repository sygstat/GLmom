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
