\name{makegrid}
\alias{makegrid}
\title{ make regular grid with square cells and round numbers }
\description{ make regular grid with square cells and round numbers }
\usage{
makegrid(x, y, n=10000, nsig=2, margin=1.05, cell.size)
}
\arguments{
\item{x}{ x-coordinate }
\item{y}{ y-coordinate }
\item{n}{ approximate number of cells in grid }
\item{nsig}{number of significant digits to which cell size and origin are rounded}
\item{margin}{margin around the x and y coordinate limits}
\item{cell.size}{cell size; if missing, a reasonable, and rounded, estimate is made}
}
\value{
data frame with the following elements:
\item{x}{ x-coordinates }
\item{y}{ y-coordinate } 
}
\examples{
data(meuse)
grd <- makegrid(meuse$x, meuse$y, 1000)
diff(grd$x)
diff(grd$y)
summary(grd)
grd <- makegrid(meuse$x, meuse$y, cell.size = 40)
diff(grd$x)
diff(grd$y)
summary(grd)
}
\keyword{dplot}
