\name{meuse.alt}
\alias{meuse.alt}
\title{Meuse river altitude data set}
\description{
This data set gives a point set with altitudes, digitized from the
1:10,000 topographical map of the Netherlands.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{x}{a numeric vector; x-coordinate (m) in RDM (Dutch topographical map coordinates) }
    \item{y}{a numeric vector; y-coordinate (m) in RDM (Dutch topographical map coordinates)}
	\item{alt}{altitude in m. above NAP (Dutch zero for sea level)}
  }
}
\usage{
data(meuse.alt)
}
\references{ \url{http:/www.gstat.org/} }
\seealso{\link{meuse}}
\keyword{datasets}
\examples{
data(meuse.alt)
xyplot(y~x, meuse.alt)
}
