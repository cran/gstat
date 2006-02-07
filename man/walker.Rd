% $Id: walker.Rd,v 1.2 2006-02-10 19:03:27 edzer Exp $
\name{walker}
\alias{walker}
\title{Walker Lake sample data set}
\description{
This is the Walker Lake sample data set (not the exhaustive
data set!), used in Isaaks and Srivastava's Applied Geostatistics.
}
\format{
  This data frame contains the following columns:
  \describe{
  \item{Id}{Identification Number}
  \item{X}{Xlocation in meter}
  \item{Y}{Ylocation in meter}
  \item{V}{V variable, concentration in ppm}
  \item{U}{U variable, concentration in ppm}
  \item{T}{T variable, indicator variable}
  }
}
\usage{
data(walker)
}
\references{ 
Applied Geostatistics
by Edward H. Isaaks, R. Mohan Srivastava;
Oxford University Press.
}
\note{
This data set was obtained from
\url{http://www.ai-geostats.org/resources/data/walker.dat};
The full (exhaustive) Walker Lake set is available from
\url{http://www.ai-geostats.org/resources/data/WalkerLake.zip}
}
\keyword{datasets}
\examples{
data(walker)
summary(walker)
}
