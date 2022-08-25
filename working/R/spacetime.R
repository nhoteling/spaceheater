# spacetime.R

spacetimer <- function(pts, dt, granularity) {
  
}


snap_to_line <- function(pts, lines) {
  pts_new <- sf::st_as_sf(maptools::snapPointsToLines(sf::as_Spatial(pts),
                                               sf::as_Spatial(lines))) 
  return(pts_new$geometry)
}

spt_add_static_interval <- function(df, interval) {
  
}

spt_interpolate <- function(df, datetime,by="day") {
  tmin <- min(as.POSIXct(datetime))
  tmax <- max(as.POSIXct(datetime))
  dt <- seq(tmin,tmax,by=by)
  
  crd <- st_coordinates(df$geometry)
  
  df.tst <- as.data.frame(approx(x=crd[,"X"],y=crd[,"Y"],n=length(dt)))
}

spt_segmentize <- function(df) {
  
}
