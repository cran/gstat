\name{fit.variogram}
\alias{fit.variogram}
\title{ Fit a Variogram Model to a Sample Variogram }
\description{
Fit a simple or nested variogram model to a sample variogram, allowing
partial fixing of parameters
}
\usage{
fit.variogram(object, model, fit.sills = T, fit.ranges = T,
	fit.method = 7, print.SSE = FALSE, debug.level = 1)
}
\arguments{
\item{object}{ sample variogram, output of \link{variogram} }
\item{model}{ variogram model, output of \link{vgm} }
\item{fit.sills}{ logical; determines whether the partial sill
coefficients (including nugget variance) should be fitted; or logical
vector: determines for each partial sill parameter whether it should be
fitted or fixed.  }
\item{fit.ranges}{ logical; determines whether the range coefficients
(excluding that of the nugget component) should be fitted; or logical
vector: determines for each range parameter whether it should be fitted
or fixed.  }
\item{fit.method}{ fitting method, used by gstat. The default method uses
weights $N_h/h^2$ with $N_h$ the number of point pairs and $h$ the
distance. This criterion is not supported by theory, but by practice. 
For other values of \code{fit.method}, see table 4.2 in the gstat
manual. }
\item{print.SSE}{ logical; if TRUE, print the (weighted) sum of squared
errors of the fitted model }
\item{debug.level}{ integer; set gstat internal debug level }
}
\value{
returns a fitted variogram model (of class \code{variogram.model})
}
\references{ \url{http://www.gstat.org/} }
\author{ Edzer J. Pebesma }
\seealso{
\link{variogram},
\link{vgm}
}
\examples{
data(meuse)
vgm1 <- variogram(log(zinc)~1, ~x+y, meuse)
fit.variogram(vgm1, vgm(1,"Sph",300,1))
}

\keyword{models}
