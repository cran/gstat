\name{plot.point.pairs}
\alias{plot.point.pairs}
\title{
Plot a point pairs, identified from a variogram cloud
}
\description{
Plot a point pairs, identified from a variogram cloud
}
\usage{
plot.point.pairs(x, data, xcol = data$x, ycol = data$y, xlab = "x coordinate", 
ylab = "y coordinate", ...) 
}
\arguments{
\item{x}{ object of class "point.pairs", obtained from the function 
\link{plot.variogram.cloud}, containing point pair indices }
\item{data}{ data frame to which the indices refer (from which the
variogram cloud was calculated) }
\item{xcol}{ numeric vector with x-coordinates of data }
\item{ycol}{ numeric vector with y-coordinates of data }
\item{xlab}{ x-axis label }
\item{ylab}{ y-axis label }
\item{...}{ arguments passed to \code{xyplot}}
}
\value{
plots the data locations, with lines connecting the point pairs identified
(and refered to by indices in) x
}
\references{ \url{http://www.gstat.org}}
\author{ Edzer J. Pebesma }
\seealso{
\link{plot.variogram.cloud}
}
\examples{
#data(meuse)
#vgm1 <- variogram(log(zinc)~1, ~x+y, meuse, cloud = TRUE)
#pp <- plot(vgm1, id = TRUE)
# identify the point pairs
#plot(pp, data = meuse) # meuse has x and y as coordinates
}

\keyword{dplot}
