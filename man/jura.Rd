% $Id: jura.Rd,v 1.2 2006-02-10 19:03:27 edzer Exp $
\name{jura}
\alias{jura}
\alias{prediction.dat}
\alias{validation.dat}
\alias{transect.dat}
\title{Jura data set}
\description{
This is the data set that comes with the Goovaerts book on
Geostatistics; it contains three data.frames:
prediction.dat, validation.dat and transect.dat.
}
\format{
  This data frame contains the following columns:
  \describe{
	\item{Xloc}{ see book }
	\item{Yloc}{ see book }
	\item{Landuse}{ see book }
	\item{Rock}{ see book }
	\item{Cd}{ see book }
	\item{Co}{ see book }
	\item{Cr}{ see book }
	\item{Cu}{ see book }
	\item{Ni}{ see book }
	\item{Pb}{ see book }
	\item{Zn}{ see book }
  }
}
\usage{
data(jura)
}
\references{ 
Goovaerts, P. 1997. Geostatistics for Natural Resources Evaluation. Oxford
Univ. Press, New-York, 483 p. Appendix C describes the Jura data set.
}
\note{
This data set was obtained from
\url{ http://home.comcast.net/~goovaerts/book.html }
}
\keyword{datasets}
\examples{
data(jura)
summary(prediction.dat)
summary(validation.dat)
summary(transect.dat)
}
