\name{meuse}
\alias{meuse}
\title{Meuse river data set}
\description{
This data set gives locations and top soil heavy metal concentrations
(ppm), along with a number of soil and landscape variables, collected in
a flood plain of the river Meuse, near the village Stein. Heavy metal
concentrations are bulk sampled from an area of approximately 15 m
x 15 m.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{Sample}{original sample number}
    \item{x}{a numeric vector; x-coordinate (m)
	in RDM (Dutch topographical map coordinates) }
    \item{y}{a numeric vector; y-coordinate (m)
	in RDM (Dutch topographical map coordinates)}
	\item{cadmium}{topsoil cadmium concentration, ppm.; note that zero cadmium
	values in the original data set have been shifted to 0.2 (half the lowest 
	non-zero value) }
	\item{copper}{topsoil copper concentration, ppm. }
	\item{lead}{topsoil lead concentration, ppm. }
	\item{zinc}{topsoil zinc concentration, ppm. } 
	\item{elev}{relative elevation}
	\item{dist}{distance to river Meuse (m)}
	\item{om}{organic matter, as percentage }
	\item{ffreq}{flooding frequency class}
	\item{soil}{soil type}
	\item{lime}{lime class}
	\item{landuse}{landuse class}
  }
}
\usage{
data(meuse)
}
\authors{ The actual field data were collected by Ruud van Rijn and
Mathieu Rikken }
\references{ 
P.A. Burrough, R.A. McDonnell, 1998. Principles of Geographical Information
Systems. Oxford University Press.

\url{http:/www.gstat.org/}
}
\author{ }
\keyword{datasets}

\examples{
data(meuse)
summary(meuse)
}
