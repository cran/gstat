\name{variogram.line}
\alias{variogram.line}
\title{ Semivariance Values For a Given Variogram Model }
\description{ Generates a semivariance values given a variogram model }
\synopsis{
variogram.line(object, maxdist, n=200, min=1.0e-6*maxdist, ids=c(0,0), dir = c(1,0,0), ...)
}
\usage{
variogram.line(object, maxdist, n, min, ids, dir, ...)
}
\arguments{
\item{object}{ variogram model for which we want semivariance function values }
\item{maxdist}{ maximum distance for which we want semivariances }
\item{n}{ number of points }
\item{min}{ minimum distance; a value slightly larger than zero is usually
used to avoid the discontinuity at distance zero if a nugget component is
present }
\item{ids}{ pair of data ids }
\item{dir}{ direction vector: unit length vector pointing the direction in 
x (East-West), y (North-South) and z (Up-Down) }
\item{...}{ignored}
}
\value{
a data frame of dimension (\code{n} x 2), with columns distance and gamma
}
\note{this function is used to plot a variogram model}
\author{ Edzer J. Pebesma }
\seealso{ \link{plot.variogram}}

\examples{
variogram.line(vgm(5, "Exp", 10, 5), 10, 10)
# anisotropic variogram, plotted in E-W direction:
variogram.line(vgm(1, "Sph", 10, anis=c(0,0.5)), 10, 10)
# anisotropic variogram, plotted in N-S direction:
variogram.line(vgm(1, "Sph", 10, anis=c(0,0.5)), 10, 10, dir=c(0,1,0))
}

\keyword{models}
