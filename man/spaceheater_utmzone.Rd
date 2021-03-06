% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spaceutils.R
\name{spaceheater_utmzone}
\alias{spaceheater_utmzone}
\title{Get UTM Zone}
\usage{
spaceheater_utmzone(lon)
}
\arguments{
\item{lon}{longitude, in GPS coordinates}
}
\value{
integer value of UTM zone
}
\description{
Get UTM Zone
}
\examples{
lon <- -73.023
spaceheater_utmzone(lon)
}
