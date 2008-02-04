% $Id: hscat.Rd,v 1.2 2008-01-07 13:25:22 edzer Exp $
\name{hscat}
\alias{hscat}
\title{
Produce h-scatterplot
}
\description{
Produces h-scatterplots, where point pairs having specific separation distances
are plotted. This function is a wrapper around xyplot.
}
\usage{
hscat(formula, data, breaks, ...)
}
\arguments{
\item{formula}{ specifies the dependent variable }
\item{data}{ data where the variable in formula is resolved }
\item{breaks}{ distance class boundaries }
\item{...}{ plotting parameters, passed to xyplot }
}
\value{
an object of class trellis; normally the h scatter plot
}
\author{ Edzer J. Pebesma }
\note{
Data pairs are plotted once, so the h-scatterplot are not symmetric.
}

\references{ \url{http://www.gstat.org/}

Pebesma, E.J., 2004. Multivariable geostatistics in S: the gstat package.
Computers \& Geosciences, 30: 683-691.

}
\examples{
data(meuse)
coordinates(meuse) = ~x+y
hscat(log(zinc)~1, meuse, c(0, 80, 120, 250, 500, 1000))
}

\keyword{models}
