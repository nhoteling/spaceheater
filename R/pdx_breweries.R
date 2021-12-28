#' Breweries near Portland, OR
#'
#' Dataset containing names and locations of breweries within 150 km of 
#' Portland, OR.
#' 
#' @format A data frame woth 251 rows and 6 variables:
#' \describe{
#'  \item{byname}{brewery name}
#'  \item{lon}{longitude}
#'  \item{lat}{latitude}
#'  \item{type}{brewery type}
#'  \item{geometry}{coordinates}
#'  \item{dist_Portland}{distance from Portland, OR, in km}
#' }
#'
#' @source derived from untappd.com in 2021, geocoded with a variety of tools
"pdx_breweries"