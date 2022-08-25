# For data handling
library(dplyr)
# library(stringr)

# For spatial things
library(sf)
library(alphahull)            # detailed pgons from point data
library(dbscan)               # spatial clustering




# Helper Functions ########################################################

#
# Density-based anomalies
# - Uses functions from spatstat
#

# Convert sf points to a density grid
points_to_density <- function(pts, sigma, buffer_size=0.15) {
  
  # Add buffer area so contours are closed
  b_check <- st_bbox(pts)
  x_buf <- buffer_size*(b_check$xmax - b_check$xmin)
  y_buf <- buffer_size*(b_check$ymax - b_check$ymin)
  buf <- ifelse(x_buf > y_buf, x_buf, y_buf)
  
  # Convert sf points to spatstat point-pattern objects
  b <- st_bbox( st_buffer(pts, buf) )
  w <- spatstat.geom::owin(xrange=c(b$xmin,b$xmax), yrange=c(b$ymin,b$ymax))
  crd <- as.data.frame(st_coordinates(pts))
  crd1 <- crd[ !duplicated(crd), ]
  x <- spatstat.geom::ppp(crd1$X, crd1$Y, window=w)
  
  # Calculate a density grid
  k <- spatstat.core::density.ppp(x, kernal="epanechnikov", sigma=sigma)
  
  return(k)
}

# TEST function: get density contours from sf points
# See outliers_by_density() function below
get_density_contours <- function(pts, sigma=0.5) {
  
  k <- points_to_density(pts, sigma)
  
  # Set up variables for contour line calculation
  xvals <- seq(from=k$xrange[1], to=k$xrange[2]-k$xstep, by=k$xstep)
  yvals <- seq(from=k$yrange[1], to=k$yrange[2]-k$ystep, by=k$ystep)
  clines <- grDevices::contourLines(xvals, yvals, t(k$v))
  
  # Convert contour lines into spatial feature lines
  d <- lapply(seq_len(length(clines)), function(i) {
    pt <- st_multipoint(matrix(c(clines[[i]]$x,clines[[i]]$y), ncol=2)) 
    #ln <- st_linestring(pt) #%>% st_sfc(crs = st_crs(pts))
    return(pt)
  })
  ln <- st_multilinestring(d) %>% st_sfc(crs = st_crs(pts))
  
  return(ln)
}

# normalize the tension to get appropriate alpha
# value for alphahull
normalize_tension <- function(pts, tension) {
  bbox <- st_bbox(pts)
  min <- st_point(c(bbox["xmin"],bbox["ymin"])) %>% st_sfc(crs=st_crs(pts))
  max <- st_point(c(bbox["xmax"],bbox["ymax"])) %>% st_sfc(crs=st_crs(pts))
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
#' Title
#'
#' @param pts 
#' @param tension 
#' @param buffer_size 
#'
#' @return
#' @export
#'
#' @examples
get_alphahull_polygon <- function(pts, tension, buffer_size=0) {
  alpha <- normalize_tension(pts, tension)
  crd <- as.data.frame(st_coordinates(pts))
  crd1 <- crd[!duplicated(crd),]
  x <- alphahull::ashape(crd1, alpha=alpha)
  
  # Build an sf polygon from alpha hull output
  df <- as.data.frame(x$edges)
  d <- lapply(seq_len(nrow(df)), function(i) {
    str_mtx <- rbind( c(df$x1[i], df$y1[i]), c(df$x2[i], df$y2[i]) )
    return( st_linestring(str_mtx) )
  })
  ln <- st_multilinestring(d)
  pgon <- st_buffer(st_polygonize(ln), buffer_size) %>%
    st_sfc(crs = st_crs(pts))
  
  # Return variable is an sf GeometryCollection
  return(pgon)
}

# Edge correction
# sometimes polygons overlap with pts that were determined
# to be 'outliers' from the algorithm.  Here, we double-check 
# for overlaps.  But, st_intersects() wants to be in planar 
# coordinates so we first transform.
edge_correction <- function(pts, pgon) {
  bbox <- st_bbox(pgon)
  crd <- c( c(mean(bbox["xmin"],bbox["xmax"])), c(mean(bbox["ymin"],bbox["ymax"])) )
  epsg_code <- spaceheater_epsg(crd[1],crd[2])
  pts1 <- pts %>% st_transform(crs=epsg_code)
  pgon1 <- pgon %>% st_transform(crs=epsg_code)
  v <- st_intersects(pts1, pgon1)
  b <- unlist(lapply(seq_len(length(v)), function(i) {x <- ifelse(length(v[[i]])>0,1,-1)}))
  return(b)
}


# Outlier Functions ###########################################################


# BOX Method
outliers_by_box <- function(pts, pct=NA) {
  # Variables -n- stuff
  #original_proj <- st_crs(pts)$input                  # get original projection
  #pts <- st_transform(pts, crs=st_crs("+proj=utm"))   # convert to UTM so we work in meters
  
  # Variables
  pct <- ifelse(!is.na(pct), pct, 0.90)
  m <- data.frame(st_coordinates(pts))
  mdx <- median(m$X)
  mdy <- median(m$Y)
  
  cntr <- st_sfc( st_point( c(mdx, mdy)), crs = st_crs(pts))
  dst <- st_distance(pts, cntr )
  
  # Spatial operations
  x <- quantile(abs(m[,1] - mdx), pct )
  y <- quantile(abs(m[,2] - mdy), pct )
  
  d.pgon <- list( rbind( c(mdx-x, mdy+y), 
                         c(mdx+x, mdy+y), 
                         c(mdx+x, mdy-y), 
                         c(mdx-x, mdy-y), 
                         c(mdx-x, mdy+y)))
  
  pgon <- st_polygon( d.pgon ) 
  
  b <- as.integer(st_within(pts, pgon))
  b[ is.na(b) ] <- -1
  
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
  m <- data.frame(st_coordinates(pts))
  md.x <- median(m$X)
  md.y <- median(m$Y)
  x <- quantile(m$X, c(lo,hi))
  y <- quantile(m$Y, c(lo,hi))
  
  # Spatial operations
  d.pgon <- list( rbind( c(x[1], y[2]), 
                         c(x[2], y[2]), 
                         c(x[2], y[1]), 
                         c(x[1], y[1]), 
                         c(x[1], y[2]) ))
  pgon <- st_polygon( d.pgon ) 
  cntr <- st_point( c(md.x, md.y) )
  b <- as.integer(st_within(pts, pgon))
  b[ is.na(b) ] <- -1
  
  # Return a list with 
  # [1] vector indicating point is inside (1) / outside (NA) of area
  # [2] polygon of area
  return(list(b,pgon))
}


# CIRCLE Method
# pts: an st_sfc( st_points ) object
# center_point: the point to measure from, if not automatically determined; should be st_point object
# r: radius, if not automatically determined
# pct: percentile from which to derive r, if not explicitly defined
#
outliers_by_circle <- function(pts, center_point=NA, r=NA, pct=NA) {
  # Variables -n- stuff
  original_proj <- st_crs(pts)$input                  # get original projection
  pts <- st_transform(pts, crs=st_crs("+proj=utm"))   # convert to UTM so we work in meters
  
  # logic to handle input variables
  if (!is.na(center_point) & !is.na(r)) {
    cntr <- st_sfc( center_point ) %>% st_set_crs( st_crs(pts)$input ) # make sure objects have the same projection
  } else {
    q <- ifelse( !is.na(pct), pct, 0.90)
    m <- data.frame(st_coordinates(pts))
    md.x <- median(m$X)
    md.y <- median(m$Y)
    cntr <- st_sfc(st_point( c(md.x, md.y))) %>% st_set_crs( st_crs(pts)$input )
  }
  
  # Spatial operations
  dst <- as.integer(st_distance(cntr, pts)) 
  r <- ifelse(!is.na(r), r, quantile(dst, c(q)))
  circle <- st_buffer(cntr, dist = r)
  b <- as.integer(st_within(pts, circle))
  b[ is.na(b) ] <- -1
  circle <- st_transform(circle, crs=st_crs(original_proj))   # convert back to orig projection
  
  # Return a list with 
  # [1] vector indicating point is inside (1) / outside (NA) of area
  # [2] polygon (circle) of area
  return(list(b,circle))
}



# POLYGON Method 
#
outliers_by_polygon <- function(pts, center_point=NA, pct=NA, tension=0.3, buffer_size=0.05) {
  
  # logic to handle input variables
  if (!is.na(center_point)) {
    cntr <- st_sfc( center_point ) %>% st_set_crs( st_crs(pts)$input ) # make sure objects have the same projection
  } else {
    q <- ifelse( !is.na(pct), pct, 0.90)
    m <- data.frame(st_coordinates(pts))
    md.x <- median(m$X)
    md.y <- median(m$Y)
    cntr <- st_sfc(st_point( c(md.x, md.y))) %>% st_set_crs( st_crs(pts)$input )
  }
  
  # Spatial operations
  dst <- as.integer(st_distance(cntr, pts))               # Compute distance from center to each point
  r <- quantile(dst, c(q))                                # distance of the specified percentile
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



# CLUSTER Method
# Note: run kNNdistplot(m, k=5) first to get a better estimate for eps
outliers_by_cluster <- function(pts, eps=0.25, MinPts=5, tension=0.5, buffer_size=0.05) {
  
  # Spatial clustering
  m <- matrix(st_coordinates(pts),ncol=2)           # Create matrix from points
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
    
    pg2 <- st_cast(pg[[1]], "POLYGON")
    return(pg2)
    #return(pg)
  })
  pgon <- st_multipolygon(d) %>% st_sfc(crs=st_crs(pts))
  
  # Adjust for polygon edge effects
  #v <- st_intersects(pts, pgon)
  #b <- unlist(lapply(seq_len(length(v)), function(i) {x <- ifelse(length(v[[i]])>0,1,-1)}))
  b <- edge_correction(pts, pgon)
  
  # Return list
  return(list(b,pgon))
}


# Check if first and last point are the same
# if not, then add a row so that first and last
# are the same, so we can convert to a polygon
# via st_polygonize()... otherwise st_polygonize
# ignores it.
#contour_to_polymatrix <- function(cln) {
#  m <- matrix(c(cln$x, cln$y), ncol=2)
#  x <- c(cln$x[1], cln$x[ length(cln$x) ])
#  y <- c(cln$y[1], cln$y[ length(cln$y) ])
#  if ( (x[1] != x[2]) | (y[1] != y[2]) ) {
#    m <- rbind(m, m[nrow(m),])
#  }
#  return(m)
#}

# DENSITY Method
# 'pts' are sf points
# 'thresh' defines the fraction of max density to use as a cutoff
# 'sigma' is passed to spatstat to determine kernel size for density calc
outliers_by_density <- function(pts, thresh=0.1, sigma=0.5) {
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
    pt <- st_multipoint(matrix(c(clines[[i]]$x,clines[[i]]$y), ncol=2)) 
    return(pt)
  })
  ln <- st_multilinestring(d) %>% st_sfc(crs = st_crs(pts))
  pgon <- st_polygonize(ln) %>% st_sfc(crs = st_crs(pts))
  
  # in case contour lines fail to close...
  if (st_is_empty(pgon)) {
    pgon <- get_alphahull_polygon(pts, tension=1.0, buffer_size = 1.5)
  }
  
  # Adjust for polygon edge effects
  #v <- st_intersects(pts, pgon)
  #b <- unlist(lapply(seq_len(length(v)), function(i) {x <- ifelse(length(v[[i]])>0,1,-1)}))
  b <- edge_correction(pts, pgon)
  
  # Return list
  return(list(b,pgon))
}


# IFOREST Method
# 'thresh' defines the quantile to use for anomaly threshold
# 
outliers_by_iforest <- function(pts, thresh=0.95, tension=0.3, buffer_size=0.05) {
  df.crd <- as.data.frame(st_coordinates(pts))
  iso = isolationForest$new(sample_size = nrow(df.crd))
  iso$fit(df.crd)
  df.scores = iso$predict(df.crd)
  q <- quantile(df.scores$anomaly_score, thresh)
  b <- ifelse(df.scores$anomaly_score>q,-1,1)
  pgon <- get_alphahull_polygon(pts[ b==1 ], tension=tension, buffer_size=buffer_size) %>%
    st_set_crs(st_crs(pts)$input)
  
  # Adjust for polygon edge effects
  #v <- st_intersects(pts, pgon)
  #b <- unlist(lapply(seq_len(length(v)), function(i) {x <- ifelse(length(v[[i]])>0,1,-1)}))
  b <- edge_correction(pts, pgon)
  
  return(list(b,pgon))
}



testjunk <- function() {
  n_bins <- 50
  mx <- max(density(df.scores$anomaly_score, n=n_bins)$y)
  p.tst <- ggplot(df.scores) +
    geom_histogram(aes(x=anomaly_score, y=..density..), bins=n_bins) +
    geom_segment(aes(x=q, xend=q,y=0,yend=0.5*mx), color="red") +
    annotate("text", x=q, y=0.5*mx, label=paste("Threshold =",round(q, digits=2)),
             hjust=-0.1, color="red") +
    theme_minimal()
  
  df.crd$val <- df.scores$anomaly_score
  df.crd$b <- ifelse(df.crd$val>q,-1,1)
  cols <- c("-1" = "#de2d26", "1" = "#009933")
  p.tst2 <- ggplot(df.crd) +
    geom_point(aes(x=X, y=Y, color=as.factor(b))) +
    scale_color_manual(values = cols) +
    theme_minimal()
}

#
# GENERIC FUNCTION CALL -- Use this one
#
# A general function to call any of the outlier methods
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