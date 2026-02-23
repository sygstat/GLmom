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
#' @source UK National River Flow Archive, Peak Flow Dataset
#'   (\url{https://nrfa.ceh.ac.uk/data/peak-flow-dataset}).
#'
#' @references 
#' Grego, J. M., Yates, P. A., & Mai, K. (2015). Standard error estimation
#' for mixed flood distributions with historic maxima. Environmetrics, 26(3), 229-242.
#' \doi{10.1002/env.2333}
#' 
#' Shin, Y., Shin, Y., & Park, J. S. (2025). Building nonstationary extreme value
#' model using L-moments. Journal of the Korean Statistical Society, 1-24.
#' \doi{10.1007/s42952-025-00325-3}
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
"PhliuAgromet"

#' Trehafod River Flow Data
#'
#' Annual maximum river flow data from the Trehafod gauging station in Wales, UK.
#' This dataset is commonly used for demonstrating non-stationary extreme value
#' analysis methods.
#'
#' @format A data frame with 53 rows and 2 columns:
#' \describe{
#'   \item{Year}{Year of observation (integer, 1968-2021)}
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
#' \donttest{
#' # Fit non-stationary GEV11 model
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
