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


#' Get EPSG Code
#' 
#' Get EPSG code for local UTM projection
#' 
#' This is a generic method to determine the EPSG code for a set of coordinates
#' \code{epsg.numeric} will take a lonlat vector whereas \code{epsg.sfc} takes an
#' sfc object and determines the EPSG code from the center of the bounding box.
#' Computation is based on a method outlined in the text Geocomputation in R,
#' by Robin Lovelace, section 7.3
#'
#' @param x Coordinates or spatial object
#'
#' @return An integer with EPSG value
#' @export
#'
#' @examples
#' lonlat <- c(-76.2, 33.1)
#' pt <- sf::st_point(crd)
#' 
#' get_epsg(crd)
#' get_epsg(pt)
get_epsg <- function(x) {UseMethod("epsg")}

# Methods
epsg.numeric <- function(x) {return(epsg1(x))}
epsg.sfc <- function(x) {
  z <- sf::st_bbox(x)
  return(epsg1( c(mean(z["xmin"],z["xmax"]),
                  mean(z["ymin"],z["ymax"]))))
}

# Calculation
epsg1 <- function(lonlat) {
  utm = (floor((lonlat[1] + 180) / 6) %% 60) + 1
  if(lonlat[2] > 0) {
    utm + 32600
  } else{
    utm + 32700
  }
}


#
# From an older gis.stackexchange post
#
epsg_old <- function(lonlat) {
  return( 32700 - round((45+lonlat[2])/90)*100 + round((183+lonlat[1])/6) )
  }


#############


#crd_to_lines
#crd_to_pgons
#crd_to_points
crd_to_points <- function(lon,lat, CRS=EPSG_WGS84) {
  d <- lapply(1:length(lon), function(i) {sf::st_point(c(lon[i],lat[i]))})
  return(d %>% sf::st_sfc(crs=CRS))
}

#
# TODO... implement these functions
#
#make_point <- function(lon, lat) { sf::st_point(c(lon, lat)) }
#make_line <- function(lonlat) { lonlat %>% as.matrix() %>% sf::st_linestring() }
#
#crd_to_segments <- function(pts, grp, thresh) {
#  if (length(pts) != length(v)) {stop("pts and v must be the same length!")}
#  data.frame(geometry=pts, v=v) %>%
#    mutate(z = purrr::map_dbl(v, function(x) {ifelse(x>=thresh, 1, 0)}),
#           g = cumsum(z)) %>%
#    group_by(g) %>%
#    summarise(m = sf::st_coordinates(geometry)) %>% tidyr::nest() %>%          
#    mutate(ln = purrr::map(data, make_line) %>% sf::st_sfc()) %>% 
#    dplyr::select(g,ln)                                                      
#}

