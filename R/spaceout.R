# For data handling
#library(dplyr)
# library(stringr)

# For spatial things
#library(sf)
#library(alphahull)            # detailed pgons from point data
#library(dbscan)               # spatial clustering




# Helper Functions ########################################################

#
# Density-based anomalies
# - Uses functions from spatstat
# - TODO: adjust for new get_epsg
points_to_density <- function(pts, sigma, buffer_size=0.15) {

  original_proj <- sf::st_crs(pts)$input                  # get original projection
  epsg_code <- get_epsg(c(mean(sf::st_coordinates(pts)[,1]),
                                mean(sf::st_coordinates(pts)[,2])))
  pts <- sf::st_transform(pts, epsg_code)   # convert to UTM so we work in meters

  # Add buffer area so contours are closed
  b_check <- sf::st_bbox(pts)
  x_buf <- buffer_size*(b_check$xmax - b_check$xmin)
  y_buf <- buffer_size*(b_check$ymax - b_check$ymin)
  buf <- ifelse(x_buf > y_buf, x_buf, y_buf)

  # Convert sf points to spatstat point-pattern objects
  b <- sf::st_buffer(sf::st_as_sfc(b_check), buf) %>% sf::st_transform(original_proj) %>% sf::st_bbox()
  w <- spatstat.geom::owin(xrange=c(b$xmin,b$xmax), yrange=c(b$ymin,b$ymax))
  crd <- as.data.frame(sf::st_coordinates(pts %>% sf::st_transform(original_proj)))
  crd1 <- crd[ !duplicated(crd), ]
  x <- spatstat.geom::ppp(crd1$X, crd1$Y, window=w)

  # Calculate a density grid
  k <- spatstat.core::density.ppp(x, kernal="epanechnikov", sigma=sigma)

  return(k)
}



# normalize the tension to get appropriate alpha
# value for alphahull
normalize_tension <- function(pts, tension) {
  bbox <- sf::st_bbox(pts)
  min <- sf::st_point(c(bbox["xmin"],bbox["ymin"])) %>% sf::st_sfc(crs=sf::st_crs(pts))
  max <- sf::st_point(c(bbox["xmax"],bbox["ymax"])) %>% sf::st_sfc(crs=sf::st_crs(pts))
  dst <- sqrt( (bbox["ymax"]-bbox["ymin"])^2 + (bbox["xmax"]-bbox["xmin"])^2)
  alpha <- as.vector(tension*dst)
  return(alpha)
}

#
# Get the alpha hull shape & convert to sf object
# 'pts' is an sf object containing points
# 'tension' controls how tightly to wrap around data points
# 'buffer_size' adds buffer around points & rounded corners
#
#' Get alpha hull polygon from collection of points
#'
#' A convenience wrapper around the alphahull::ashape() function.
#' The function creates a polygon around a collection of sf points, given inputs
#' for tension and an optional buffer.
#'
#' @param pts sf points
#' @param tension value between 0 and 1, controls how tight polygon should wrap around points;
#' value is fraction of diagonal of bounding box
#' @param buffer_size add an optional buffer around points; value in coordinate units
#'
#' @return sf GeometryCollection object with a polygon
#' @export
#'
#' @seealso \code{\link{get_alphahull_polygon()}}
#'
#' @examples
#' data("pdx_breweries")
#' pts <- pdx_breweries$geometry
#' x1 <- get_alphahull_polygon(pts)
#' x2 <- get_alphahull_polygon(pts, tension=0.2, buffer_size=0.1)
#' plot(pts)
#' plot(x1, add=TRUE)
#' plot(x2, add=TRUE)
get_alphahull_polygon <- function(pts, tension=0.5, buffer_size=0) {
  alpha <- normalize_tension(pts, tension)
  crd <- as.data.frame(sf::st_coordinates(pts))
  crd1 <- crd[!duplicated(crd),]
  x <- alphahull::ashape(crd1, alpha=alpha)

  # Build an sf polygon from alpha hull output
  df <- as.data.frame(x$edges)
  d <- lapply(seq_len(nrow(df)), function(i) {
    str_mtx <- rbind( c(df$x1[i], df$y1[i]), c(df$x2[i], df$y2[i]) )
    return( sf::st_linestring(str_mtx) )
  })
  ln <- sf::st_multilinestring(d)
  pgon <- sf::st_buffer(sf::st_polygonize(ln), buffer_size) %>%
    sf::st_sfc(crs = sf::st_crs(pts))

  # Return variable is an sf GeometryCollection
  return(pgon)
}

# Edge correction
# sometimes polygons overlap with pts that were determined
# to be 'outliers' from the algorithm.  Here, we double-check
# for overlaps.  But, st_intersects() wants to be in planar
# coordinates so we first transform.
edge_correction <- function(pts, pgon) {
  #bbox <- sf::st_bbox(pgon)
  #crd <- c( c(mean(bbox["xmin"],bbox["xmax"])), c(mean(bbox["ymin"],bbox["ymax"])) )
  #epsg_code <- spaceheater_epsg(crd[1],crd[2])
  epsg_code <- get_epsg(pgon)
  pts1 <- pts %>% sf::st_transform(crs=epsg_code)
  pgon1 <- pgon %>% sf::st_transform(crs=epsg_code)
  v <- sf::st_intersects(pts1, pgon1)
  b <- unlist(lapply(seq_len(length(v)), function(i) {x <- ifelse(length(v[[i]])>0,1,-1)}))
  return(b)
}


# Outlier Functions ###########################################################


#
#' Spatial outliers: BOX Method
#'
#' Spatial outlier method that uses distance from a central point to draw a box
#' around a collection of points.  The distance to box boundaries is determined from the 90th
#' percentile, or by user input.
#' TODO: add manual center-point option, different xy distances
#' TODO: adjust code for new epsg function
#'
#' @param pts sf points
#' @param pct (optional) quantile to determine distance from center
#'
#' @return list with outlier status and sf polygon
#' @export
#'
#' @family outlier functions
#'
#' @examples
#' data("pdx_breweries")
#' pts <- pdx_breweries$geometry
#' x <- outliers_by_box(pts)
#' plot(pts)
#' plot(pts[x[[1]]==-1], col="red", add=TRUE)
#' plot(x[[2]], add=TRUE)
outliers_by_box <- function(pts, pct=NA) {
  # Variables -n- stuff
  original_proj <- sf::st_crs(pts)$input                  # get original projection
  epsg_code <- get_epsg(c(mean(sf::st_coordinates(pts)[,1]),
                                mean(sf::st_coordinates(pts)[,2])))
  pts <- sf::st_transform(pts, epsg_code)   # convert to UTM so we work in meters

  # Variables
  pct <- ifelse(!is.na(pct), pct, 0.90)
  m <- data.frame(sf::st_coordinates(pts))
  mdx <- stats::median(m$X)
  mdy <- stats::median(m$Y)

  cntr <- sf::st_sfc( sf::st_point( c(mdx, mdy)), crs = sf::st_crs(pts))
  dst <- sf::st_distance(pts, cntr )

  # Spatial operations
  x <- stats::quantile(abs(m[,1] - mdx), pct )
  y <- stats::quantile(abs(m[,2] - mdy), pct )

  d.pgon <- list( rbind( c(mdx-x, mdy+y),
                         c(mdx+x, mdy+y),
                         c(mdx+x, mdy-y),
                         c(mdx-x, mdy-y),
                         c(mdx-x, mdy+y)))

  pgon <- sf::st_polygon( d.pgon ) %>% sf::st_sfc(crs=epsg_code)

  b <- as.integer(sf::st_within(pts, pgon))
  b[ is.na(b) ] <- -1

  pgon <- sf::st_transform(pgon, crs=sf::st_crs(original_proj))   # convert back to orig projection
  # Return a list with
  # [1] vector indicating point is inside (1) / outside (NA) of area
  # [2] polygon of area
  return(list(b,pgon))
}


# This was the original box method...
outliers_by_quantile <- function(pts, pct_lo=NA, pct_hi=NA) {
  # Variables: Use Interquartile Range by default
  lo <- ifelse(!is.na(pct_lo), pct_lo, 0.25)
  hi <- ifelse(!is.na(pct_hi), pct_hi, 0.75)
  m <- data.frame(sf::st_coordinates(pts))
  md.x <- stats::median(m$X)
  md.y <- stats::median(m$Y)
  x <- stats::quantile(m$X, c(lo,hi))
  y <- stats::quantile(m$Y, c(lo,hi))

  # Spatial operations
  d.pgon <- list( rbind( c(x[1], y[2]),
                         c(x[2], y[2]),
                         c(x[2], y[1]),
                         c(x[1], y[1]),
                         c(x[1], y[2]) ))
  pgon <- sf::st_polygon( d.pgon )
  cntr <- sf::st_point( c(md.x, md.y) )
  b <- as.integer(sf::st_within(pts, pgon))
  b[ is.na(b) ] <- -1

  # Return a list with
  # [1] vector indicating point is inside (1) / outside (NA) of area
  # [2] polygon of area
  return(list(b,pgon))
}



#' Spatial outliers: CIRCLE Method
#'
#' Spatial outlier method that uses distance from a central point to draw a bounding circle
#' from which to determine outlier status.  Use this method to isolate points within some
#' distance r from a given point.  The method will automatically compute a center point
#' and distance if none are provided, based on an input percentile, or 0.90 if nothing is provided
#'
#' @param pts sf points
#' @param center_point (optional) sf point, the point to measure from
#' @param r (optional) radius around center point, in meters
#' @param pct (optional) percentile to use for determining r, if not provided; 0.90 by default
#'
#' @return list with outlier status and sf polygon
#' @export
#'
#' @family outlier functions
#'
#' @examples
#' data("pdx_breweries")
#' pts <- pdx_breweries$geometry
#' x <- outliers_by_circle(pts)
#' plot(pts)
#' plot(pts[x[[1]]==-1], col="red", add=TRUE)
#' plot(x[[2]], add=TRUE)
outliers_by_circle <- function(pts, center_point=NA, r=NA, pct=NA) {
  # Variables -n- stuff
  original_proj <- sf::st_crs(pts)$input                  # get original projection
  epsg_code <- get_epsg(c(mean(sf::st_coordinates(pts)[,1]),
                                mean(sf::st_coordinates(pts)[,2])))
  pts <- pts %>% sf::st_transform(epsg_code)   # convert to UTM so we work in meters

  # logic to handle input variables
  if (!is.na(center_point) & !is.na(r)) {
    cntr <- sf::st_sfc( center_point ) %>% sf::st_set_crs( sf::st_crs(pts)$input ) # make sure objects have the same projection
  } else {
    q <- ifelse( !is.na(pct), pct, 0.90)
    m <- data.frame(sf::st_coordinates(pts))
    md.x <- stats::median(m$X)
    md.y <- stats::median(m$Y)
    cntr <- sf::st_sfc(sf::st_point( c(md.x, md.y))) %>% sf::st_set_crs( sf::st_crs(pts)$input )
  }

  # Spatial operations
  dst <- as.integer(sf::st_distance(cntr, pts))
  r <- ifelse(!is.na(r), r, stats::quantile(dst, c(q)))
  circle <- sf::st_buffer(cntr, dist = r)
  b <- as.integer(sf::st_within(pts, circle))
  b[ is.na(b) ] <- -1
  circle <- sf::st_transform(circle, crs=sf::st_crs(original_proj))   # convert back to orig projection

  # Return a list with
  # [1] vector indicating point is inside (1) / outside (NA) of area
  # [2] polygon (circle) of area
  return(list(b,circle))
}



#
#
#' Spatial outliers: POLYGON Method
#'
#' Spatial outlier method that uses distance from a central point to draw a bounding polygon
#' from which to determine outlier status.  The outliers should be identical to the circle
#' method, except the shape will account for anisotropy in the the point pattern.  The polygon
#' is built with alpha hull, so the user can specify tension and buffer size.  Note that, for
#' very low (tighter) tension values edge effects may lead to differences from the circle method.
#'
#' @param pts sf points
#' @param center_point (optional) sf point with location to measure from
#' @param pct (optional) percentile to use for determining distance from center point
#' @param tension (optional) value from 0 to 1 to control how tight polygon should wrap around points
#' @param buffer_size (optional) buffer around exterior points, in coordinate units
#'
#' @return list with outlier status and sf polygon
#' @export
#'
#' @seealso \code{\link{get_alphahull_polygon}}
#' @family outlier functions
#'
#' @examples
#' data("pdx_breweries")
#' pts <- pdx_breweries$geometry
#' x <- outliers_by_polygon(pts)
#' plot(pts)
#' plot(pts[x[[1]]==-1], col="red", add=TRUE)
#' plot(x[[2]], add=TRUE)
outliers_by_polygon <- function(pts, center_point=NA, pct=NA, tension=0.3, buffer_size=0.05) {

  # logic to handle input variables
  if (!is.na(center_point)) {
    cntr <- sf::st_sfc( center_point ) %>% sf::st_set_crs( sf::st_crs(pts)$input ) # make sure objects have the same projection
  } else {
    q <- ifelse( !is.na(pct), pct, 0.90)
    m <- data.frame(sf::st_coordinates(pts))
    md.x <- stats::median(m$X)
    md.y <- stats::median(m$Y)
    cntr <- sf::st_sfc(sf::st_point( c(md.x, md.y))) %>% sf::st_set_crs( sf::st_crs(pts)$input )
  }

  # Spatial operations
  dst <- as.integer(sf::st_distance(cntr, pts))               # Compute distance from center to each point
  r <- stats::quantile(dst, c(q))                                # distance of the specified percentile
  b <- ifelse(dst <= r, +1, -1)                           # which points are within distance r from center

  # alpha hull business #######################################################
  # Get coordinates from points & remove any duplicates
  pgon <- get_alphahull_polygon(pts[ b==1 ], tension=tension, buffer_size=buffer_size)

  # Adjust for polygon edge effects
  #v <- st_intersects(pts, pgon)
  #b <- unlist(lapply(seq_len(length(v)), function(i) {x <- ifelse(length(v[[i]])>0,1,-1)}))
  b <- edge_correction(pts, pgon)

  # Return a list with ########################################################
  # [1] vector indicating point is inside (1) / outside (NA) of area
  # [2] polygon (circle) of area
  return(list(b,pgon))
}




#' Spatial outliers: CLUSTER Method
#'
#' Spatial outlier method based on the dbscan spatial clustering algorithm.  User should
#' run dbscan::kNNdistplot() first to get a reasonable value for eps.  The dbscan algorithm
#' clusters according to a nominal distance between points (eps) and minimum number of
#' points per cluster (MinPts) and puts all points that don't neatly fit into a group
#' into cluster zero.  This method draws a polygon around the non-cluster-zero points,
#' applies some edge corrections, and determines outlier status accordingly.
#'
#' @param pts sf points
#' @param eps parameter passed to dbscan, distance between neighboring points
#' @param MinPts parameter passed to dbscan, min points per cluster
#' @param tension (optional) value from 0 to 1 to control how tight polygon should wrap around points
#' @param buffer_size (optional) buffer around exterior points, in coordinate units
#'
#' @return list with outlier status and sf polygon
#' @export
#'
#' @seealso \code{\link{get_alphahull_polygon()}}
#' @family outlier functions
#'
#' @examples
#' data("pdx_breweries")
#' pts <- pdx_breweries$geometry
#' x <- outliers_by_cluster(pts)
#' plot(pts)
#' plot(pts[x[[1]]==-1], col="red", add=TRUE)
#' plot(x[[2]], add=TRUE)
outliers_by_cluster <- function(pts, eps=0.25, MinPts=5, tension=0.5, buffer_size=0.05) {

  # Spatial clustering
  m <- matrix(sf::st_coordinates(pts),ncol=2)           # Create matrix from points
  dbsc <- dbscan::dbscan(m, eps=eps, minPts = MinPts)   # Run dbscan

  # An alternative is to use hdbscan so that eps is not needed, but
  # in the examples I tested this tends to produce too much fine structure
  #hdbsc <- dbscan::hdbscan(m, minPts = MinPts)


  x <- dbsc$cluster                                 # Get clusters from dbscan

  # We can subtract the "outliers" from the polygon...
  # ...or just deal with it...
  # below, we use intersect to define keep/toss instead.
  #pg0 <- get_alphahull_polygon(pts[ x==0 ], tension=0.3)

  # Get polygons from clustering results
  d <- lapply(1:max(x), function(i) {
    pg <- get_alphahull_polygon(pts[ x==i ], tension=tension, buffer_size=buffer_size)

    # If subtracting "outlier" shape then uncomment this
    #pg1 <- st_difference(pg,pg0)
    #pg2 <- st_cast(pg1[[1]], "POLYGON")  # Not sure why this was needed

    pg2 <- sf::st_cast(pg[[1]], "POLYGON")
    return(pg2)
    #return(pg)
  })
  pgon <- sf::st_multipolygon(d) %>% sf::st_sfc(crs=sf::st_crs(pts))

  # Adjust for polygon edge effects
  #v <- st_intersects(pts, pgon)
  #b <- unlist(lapply(seq_len(length(v)), function(i) {x <- ifelse(length(v[[i]])>0,1,-1)}))
  b <- edge_correction(pts, pgon)

  # Return list
  return(list(b,pgon))
}





#' Spatial outliers: DENSITY Method
#'
#' Spatial outlier method that uses kernel density to determine outlier status.
#' The algorithm uses spatstat to compute a density grid, grDevices to create
#' contour lines based on user-defined threshold value, and then determines outliers
#' by intersecting points with the contours.  The threshold value is based on some
#' fraction of the maximum density, 0.1 in the absence of user input.
#'
#' @param pts sf points
#' @param thresh fraction of maximum density to use as outlier threshold
#' @param sigma value passed to spatstat for kernel density size;
#' the value is used as standard deviation of Gaussian for smoothing
#'
#' @return list with outlier status and sf polygon
#' @export
#'
#' @family outlier functions
#' @seealso \code{\link{points_to_density()}}
#' @examples
#' data("pdx_breweries")
#' pts <- pdx_breweries$geometry
#' x <- outliers_by_density(pts)
#' plot(pts)
#' plot(pts[x[[1]]==-1], col="red", add=TRUE)
#' plot(x[[2]], add=TRUE)
outliers_by_density <- function(pts, thresh=0.1, sigma=0.1) {
  k <- points_to_density(pts, sigma)  # Function defined above
  zlevel <- thresh*(max(k$v))
  # Set up variables for contour line calculation
  xvals <- seq(from=k$xrange[1], to=k$xrange[2]-k$xstep, by=k$xstep)
  yvals <- seq(from=k$yrange[1], to=k$yrange[2]-k$ystep, by=k$ystep)
  clines <- grDevices::contourLines(xvals, yvals, t(k$v), levels=zlevel)

  # Convert contour lines into spatial feature lines
  d <- lapply(seq_len(length(clines)), function(i) {
    #mtx <- contour_to_polymatrix(clines[[i]])
    #pt <- st_multipoint(mtx)
    pt <- sf::st_multipoint(matrix(c(clines[[i]]$x,clines[[i]]$y), ncol=2))
    return(pt)
  })
  ln <- sf::st_multilinestring(d) %>% sf::st_sfc(crs = sf::st_crs(pts))
  pgon <- sf::st_polygonize(ln) %>% sf::st_sfc(crs = sf::st_crs(pts))

  # in case contour lines fail to close...
  if (sf::st_is_empty(pgon)) {
    pgon <- get_alphahull_polygon(pts, tension=1.0, buffer_size = 1.5)
  }

  # Adjust for polygon edge effects
  #v <- st_intersects(pts, pgon)
  #b <- unlist(lapply(seq_len(length(v)), function(i) {x <- ifelse(length(v[[i]])>0,1,-1)}))
  b <- edge_correction(pts, pgon)

  # Return list
  return(list(b,pgon))
}



#' Spatial outliers: IFOREST Method
#'
#' Spatial outliers based on the Isolation Forest algorithm implemented in the
#' solitude package.  IForest is a tree-based method used to determine nominal difference
#' from normal values; it is not traditionally applied to spatial data, so this
#' method is really just experimental.  Threshold value defines the percentile to
#' use for determining outlier status.
#'
#' @param pts sf points
#' @param thresh percentile value to use for determining outlier status
#' @param tension (optional) value from 0 to 1 to control how tight polygon should wrap around points
#' @param buffer_size (optional) buffer around exterior points, in coordinate units
#'
#' @return list with outlier status and sf polygon
#' @export
#'
#' @family outlier functions
#'
#' @examples
#' data("pdx_breweries")
#' pts <- pdx_breweries$geometry
#' x <- outliers_by_iforest(pts)
#' plot(pts)
#' plot(pts[x[[1]]==-1], col="red", add=TRUE)
#' plot(x[[2]], add=TRUE)
outliers_by_iforest <- function(pts, thresh=0.95, tension=0.3, buffer_size=0.05) {
  df.crd <- as.data.frame(sf::st_coordinates(pts))
  iso = solitude::isolationForest$new(sample_size = nrow(df.crd))
  iso$fit(df.crd)
  df.scores = iso$predict(df.crd)
  q <- stats::quantile(df.scores$anomaly_score, thresh)
  b <- ifelse(df.scores$anomaly_score>q,-1,1)
  pgon <- get_alphahull_polygon(pts[ b==1 ], tension=tension, buffer_size=buffer_size) %>%
    sf::st_set_crs(sf::st_crs(pts)$input)

  # Adjust for polygon edge effects
  #v <- st_intersects(pts, pgon)
  #b <- unlist(lapply(seq_len(length(v)), function(i) {x <- ifelse(length(v[[i]])>0,1,-1)}))
  b <- edge_correction(pts, pgon)

  return(list(b,pgon))
}





#' Spatial outliers
#'
#' A convenience wrapper used to call any of the spatial outlier methods implemented
#' in \code{spaceheater}.
#'
#' @param pts sf points
#' @param method which method to use, "box", "circle", "polygon", "density", "cluster", "iforest"
#' @param ... input parameters specific to individual methods
#'
#' @return list with outlier status and sf polygon
#' @export
#'
#' @family outlier functions
#'
#' @examples
#' data("pdx_breweries")
#' pts <- pdx_breweries$geometry
#' x <- spatial_outliers(pts,method="polygon")
#' plot(pts)
#' plot(pts[x[[1]]==-1], col="red", add=TRUE)
#' plot(x[[2]], add=TRUE)
spatial_outliers <- function(pts, method = "polygon", ...) {
  d <- switch(method,
              "box" = outliers_by_box(pts, ...),
              "circle" = outliers_by_circle(pts, ...),
              "polygon" = outliers_by_polygon(pts, ...),
              "cluster" = outliers_by_cluster(pts, ...),
              "quantile" = outliers_by_quantile(pts, ...),
              "density" = outliers_by_density(pts, ...),
              "iforest" = outliers_by_iforest(pts, ...))
  # Return variable is a list containing two items:
  # 1. Vector indicating outlier (-1) or not (+1)
  # 2. Polygon outlining the relevant area
  return(d)
}
