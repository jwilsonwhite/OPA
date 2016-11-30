#' Map of semi-desert grassland plants
#' 
#' A spatial dataset of class \code{"ppp"} containing the locations and other information collected on semi-arid grassland plants 
#' in 178 permanent 1x1 m quadrats at the Santa Rita Experimental Range, AZ, USA. Data were mapped annually from 1915-1933 and again in 1947.
#' Details can be found in Anderson et al. (2012).
#' Objects of class \code{"ppp"} are essentially data frames; the variables are as follows:
#' 
#' \itemize{
#' \item \code{x}. \emph{x}-coordinate of the centroid of each plant
#' \item \code{y}. \emph{y}-coordinate of the centroid of each plant
#' \item \code{marks}. Data frame containing additional data for each plant. Columns are
#' \enumerate{
#' \item radius. The radius of each plant rosette, calculated from perimeter (assuming a circular shape)
#' \item quad. The quadrat identifier code
#' \item year. The year sampled, in YY format ('15' = 1915)
#' \item Species. The species name
#' \item Clone. The clone status of the plant
#' \item Seeding. Is the plant a seedling? (\code{Y} or \code{N})
#' \item Area. The basal area in square meters
#' \item Length. The length of the rosette perimeter, in meters
#' }}
#' 
#'
#' @docType data
#' @keywords datasets
#' @name grasslandplants
#' @usage data(grasslandplants)
#' @format a data frame with 133554 rows and 9 variables
#' @source \href{http://onlinelibrary.wiley.com/doi/10.1890/11-2200.1/full}{Anderson JA, McClaran MP, Adler PB (2012) Cover and density of semi-desert grassland plants in permanent quadrats mapped from 1915 to 1947. Ecology 93:1492} 
NULL

