% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spaceout.R
\name{outliers_by_circle}
\alias{outliers_by_circle}
\title{Spatial outliers: CIRCLE Method}
\usage{
outliers_by_circle(pts, center_point = NA, r = NA, pct = NA)
}
\arguments{
\item{pts}{sf points}

\item{center_point}{(optional) sf point, the point to measure from}

\item{r}{(optional) radius around center point, in meters}

\item{pct}{(optional) percentile to use for determining r, if not provided; 0.90 by default}
}
\value{
list with outlier status and sf polygon
}
\description{
Spatial outlier method that uses distance from a central point to draw a bounding circle
from which to determine outlier status.  Use this method to isolate points within some
distance r from a given point.  The method will automatically compute a center point
and distance if none are provided, based on an input percentile, or 0.90 if nothing is provided
}
\examples{
data("pdx_breweries")
pts <- pdx_breweries$geometry
x <- outliers_by_circle(pts)
plot(pts)
plot(pts[x[[1]]==-1], col="red", add=TRUE)
plot(x[[2]], add=TRUE)
}
\seealso{
Other outlier functions: 
\code{\link{outliers_by_box}()},
\code{\link{outliers_by_cluster}()},
\code{\link{outliers_by_density}()},
\code{\link{outliers_by_iforest}()},
\code{\link{outliers_by_polygon}()},
\code{\link{spatial_outliers}()}
}
\concept{outlier functions}
