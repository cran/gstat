\name{fulmar}
\alias{fulmar}
\title{Fulmaris glacialis data}
\description{
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
\authors{ Dutch National Institute for Coastal and Marine Management (RIKZ),
\url{http://www.rikz.nl/}}
\references{ 
}
\note{ }
\seealso{\link{ncp.grid}}
\keyword{datasets}
\examples{
data(fulmar)
summary(fulmar)
}
