\name{select.spatial}
\alias{select.spatial}
\title{ select points spatially }
\description{
select a number of points by digitizing the area they fall in
}
\usage{
select.spatial(x = data$x, y = data$y, data, pch = "+")
}
\arguments{
\item{x}{numerical array of x-coordinates of points}
\item{y}{numerical array of y-coordinates of points}
\item{data}{optional; data frame containing variables x and y}
\item{pch}{plotting character to be used for points}
}
\value{ array with indexes (row numbers) of points inside
the polygon digitized }
\seealso{\link{point.in.polygon}, \link{locator}}
\examples{
data(meuse)
## the following command requires user interaction: left mouse
## selects points, right mouse ends digitizing
# select.spatial(data=meuse)
}
\keyword{models}
