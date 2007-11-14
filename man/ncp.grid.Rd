% $Id: ncp.grid.Rd,v 1.4 2007-11-13 21:34:26 edzer Exp $
\name{ncp.grid}
\alias{ncp.grid}
\title{Grid for the NCP, the Dutch part of the North Sea}
\description{
Gridded data for the NCP, the Dutch part of the North Sea, for
a 5 km x 5 km grid; stored as data.frame.
}
\format{
  This data frame contains the following columns:
  \describe{
   \item{x}{x-coordinate, UTM31}
   \item{y}{y-coordinate, UTM31}
   \item{depth}{sea water depth, m}
   \item{coast}{distance to coast, m}
   \item{area}{identifier for administrative sub-areas}
  }
}
\usage{
data(ncp.grid)
}
\author{Dutch National Institute for Coastal and Marine Management (RIKZ);
data compiled for R by Edzer J. Pebesma }
\references{ 
}
\note{ }
\seealso{\link{fulmar}}
\keyword{datasets}
\examples{
data(ncp.grid)
summary(ncp.grid)
}
