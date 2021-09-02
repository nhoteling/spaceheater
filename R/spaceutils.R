# spaceheater.R
# Miscellaneous functions here
#


# Get UTM Zone & EPSG code ########################
# UTM based on:
# https://gis.stackexchange.com/questions/209267/r-return-the-utm-zone-that-a-wgs84-point-belongs-to
# EPSG based on:
# https://gis.stackexchange.com/questions/190198/how-to-get-appropriate-crs-for-a-position-specified-in-lat-lon-coordinates/190209#190209
#
# To check:
# https://mangomap.com/robertyoung/maps/69585/what-utm-zone-am-i-in-#
#
# OR
# > epsg <- rgdal::make_EPSG()
# > epsg %>% filter(code == spaceheater_epsg(lon,lat))
# The result should make sense
#

EPSG_WGS84 <- 4326

#
#' Get UTM Zone
#'
#' @param lon longitude, in GPS coordinates
#'
#' @return integer value of UTM zone
#' @export
#'
#' @examples
#' lon <- -73.023
#' spaceheater_utmzone(lon)
spaceheater_utmzone <- function(lon) {
  return( floor((lon + 180) / 6) + 1 )
  #return( round(183+lon)/6 )
}


#' Get EPSG Code for local UTM projection
#'
#' Compute the EPSG code for local UTM projection, based on GPS inputs.
#' Calculation is based on standard formula:
#' 32700 - round((45+lat)/90)*100 + round((183+lon)/6)
#'
#' @param lon longitude GPS coordinate
#' @param lat latitude GPS coordinate
#'
#' @return integer value for EPSG code
#' @export
#'
#' @examples
#' crd <- c(-76.3451, 34.3352)
#' spaceheater_epsg(crd[1], crd[2])
spaceheater_epsg <- function(lon, lat) {
  return( 32700 - round((45+lat)/90)*100 + round((183+lon)/6) )
}

#######
