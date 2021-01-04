
#' Data example for SST prediction
#'
#' A dataset containing the (not postprocessed) predictions by NorCPM and the corresponding observations (OISST) for monthly mean sea surface temperature.
#' The data is restricted to January and February over the Northern Atlantic Ocean, in order to keep it at a moderate size.
#' See the main body of the paper for more details. For visualization of the data please check the function \code{plot_diagnostic}.
#' For getting access to the full data containing all months and the entire globe please contact the authors of the article.
#' 
#' 
#'
#' @format A data frame with 557.700 rows and 17 variables:
#' \describe{
#'   \item{Lon, Lat}{Longitude and Latitude.}
#'   \item{Ens1,...,Ens9}{NorCPM prediction ensemble members. The predictions are in degree celsius, and are predictions for the monthly mean SST of the month specified by the columns month, year. Missing values indicate land. The forecasts are initialized at January first of the same year, see the paper for more details.}
#'   \item{Ens_bar, Ens_sd}{Ensemble mean and standard deviation.}
#'   \item{SST_bar}{Observed monthly mean SST in degree celsius provided by OISST.}
#'   \item{year,month}{year and month of the observation and for which the predictions are issued.}
#'   \item{grid_id}{An identifyer of the grid point.}
#' }
#'
#' @source The OISST data is available at \url{https://www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/access/.} NorCPM is currently under development and the predictions used in this work are not publicly available online, but will be made available upon request. Please contact one of the authors.
#' 
"DT"