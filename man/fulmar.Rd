% $Id: fulmar.Rd,v 1.7 2008-07-03 11:49:08 edzer Exp $
\name{fulmar}
\alias{fulmar}
\title{Fulmaris glacialis data}
\description{ Airborne counts of Fulmaris glacialis during the Aug/Sept 1998 and 1999
flights on the Netherlands part of the North Sea (NCP).
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{year}{year of measurement: 1998 or 1999}
	\item{x}{x-coordinate in UTM31} 
	\item{y}{y-coordinate in UTM31}
	\item{depth}{sea water depth, in m}
	\item{coast}{distance to coast, in m}
	\item{fulmar}{observed density (number of birds per square km)}
  }
}
\usage{
data(fulmar)
}
\author{ Dutch National Institute for Coastal and Marine Management (RIKZ),
\url{http://www.rikz.nl/}
}
\seealso{\link{ncp.grid}

E.J. Pebesma, R.N.M. Duin, P.A. Burrough, 2005. Mapping Sea Bird
Densities over the North Sea: Spatially Aggregated Estimates and Temporal
Changes. Environmetrics 16, (6), p 573-587.
}
\keyword{datasets}
\examples{
data(fulmar)
summary(fulmar)
}
