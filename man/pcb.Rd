\name{pcb}
\alias{pcb}
\title{PCB138 measurements in sediment at the NCP, the Dutch part of the North Sea}
\description{
This data set gives a point set with altitudes, digitized from the
1:10,000 topographical map of the Netherlands.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{year}{ measurement year }
    \item{x}{ x-coordinate; UTM31 }
    \item{y}{ y-coordinate; UTM31 }
    \item{coast}{ distance to coast, m.}
    \item{depth}{ sea water depth, m. }
    \item{PCB138}{ PCB-138, measured on the sediment fraction smaller than
	63 $\mu$m, in $\mu$g/kg dry matter; BUT SEE NOTE BELOW }
	\item{yf}{ year; as factor}
  }
}

\usage{
data(pcb)
}
\note{ A note of caution: The PCB-138 data are provided only to be able
to re-run the analysis done in Pebesma and Duin (2004; see references
below). If you want to use these data for comparison with PCB measurements
elsewhere, or if you want to compare them to regulation standards, or
want to use these data for any other purpose, you should first contact
\url{mailto:basisinfodesk@rikz.rws.minvenw.nl}.  The reason for this is
that several normalisations were carried out that are not reported here,
nor in the paper below. }

\references{ \url{http:/www.gstat.org/}, \url{http://www.rikz.nl/} 

Edzer J. Pebesma, Richard N.M. Duin, 2004.  Spatio-temporal mapping of sea
floor sediment pollution in the North Sea. Paper presented at GeoENV2004,
Oct 12-14, 2004, Neuchatel; proceedings to be published by Springer. A copy
of the paper can be requested from \url{mailto:e.pebesma@geog.uu.nl}
}

\seealso{\link{ncp.grid}}
\keyword{datasets}
\examples{
data(pcb)
library(lattice)
xyplot(y~x|as.factor(yf), pcb, aspect = "iso")
# demo(pcb)
}
