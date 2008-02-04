% $Id: fulmar.Rd,v 1.6 2007-11-16 12:59:47 edzer Exp $
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
\url{http://www.rikz.nl/}}
\seealso{\link{ncp.grid}}
\keyword{datasets}
\examples{
data(fulmar)
summary(fulmar)
}
